import pandas
from meandist2 import build_knn_dict, mean_dist_cohort
from sklearn.model_selection import GridSearchCV
import json
import optparse
import sys
from sklearn.base import BaseEstimator, ClassifierMixin
import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import logging
from tqdm import tqdm

from pedigree import *

logging.basicConfig(
    format='%(asctime)s %(levelname)-8s %(message)s',
    level=logging.INFO,
    datefmt='%Y-%m-%d %H:%M:%S')

'''
Generate test_case.txt from ped-aware/implementation/match_against_real_ibds.py
Sample command used in ped-aware/knn-dist-min/algos: python3 false_positive_percentage_removal_test.py -s ../extended_pedigree_final.fam -t test_case1.txt
'''

def parse_args(description):
    parser = optparse.OptionParser(description=description)
    parser.add_option("-s", "--struct_filename", type="string", \
                      help="input .txt file with pedigree structure")
    parser.add_option("-t", "--testcase_filename", type="string", \
                      help="input .txt file with test cases for cohorts")
    parser.add_option("-o", "--output_filename", type="string", \
                      help="output .txt file with pre and post-processed cohort \
                        with matching ibd length")

    mandatories = ["struct_filename", "testcase_filename", "output_filename"]
    (opts, args) = parser.parse_args()
    for m in mandatories:
        if not opts.__dict__[m]:
            parser.print_help()
            sys.exit()
    return opts

def build_dataset_for_gridsearch(row_data):
    # row_data comes from test_case.txt. It represents one row of data that contains two columns: cohort candidates (from germline or TPBWT) and
    # true cohort members. The different members in both columns are separated by semi-colons
    data = dict()
    columns = row_data.strip().split()
    software_ibd_length = columns[0]
    cohort_candidates = [x for x in columns[1].split(';')]
    true_cohort_ibd_length = columns[2]
    true_cohort_members = [x for x in columns[3].split(';')]
    return software_ibd_length, cohort_candidates, true_cohort_ibd_length, true_cohort_members

def read_test_cases(test_file):

    test_data = test_file.readlines()
    test_file.close()
    test_cases = test_data[1:]
    return test_cases


def scoring_func_for_gridsearch(func_output, desired_output):
    # this function will evaluate how closely the function output matches with the desired output.
    # the function output in this case will be the predicted true cohort members while the desired output
    # will be the actual true cohort members. This function will return the Jaccard similarity score
    # calculated_val = set(func_output)
    # expected_val = set(desired_output)
    # return len(calculated_val & expected_val) / len(calculated_val | expected_val)
    
    # Instead of using the Jaccard similarity score, let us use a metric that is focused on the number
    # of false positives since that's what the algorithm was ultimately designed for.
    calculated_val = set(func_output)
    expected_val = set(desired_output)
    false_positives = calculated_val - expected_val
    return len(false_positives) / len(calculated_val)


def calculate_fp_tp_fn(cohort, true_cohort_members):
    false_positives = len(cohort - true_cohort_members)
    true_positives = len(cohort & true_cohort_members)
    false_negatives = len(true_cohort_members - cohort)

    return false_positives, true_positives, false_negatives


if __name__ == "__main__":
    args = parse_args("pedigree args")
    test_file_prefix = args.testcase_filename.split('/')[-1].split('txt')[0][:-1]

    single_edge_dict = read_ped_file(args.struct_filename)
    logging.info("finish counting edges between all pairs of individuals")
    test_file = open(args.testcase_filename, 'r')
    test_cases = read_test_cases(test_file)
    logging.info("finish reading test cases")

    param_grid = {
        'num_neighbors': [2, 8, 14, 20, 26, 32, 38, 44, 50], 
        't': [0.00625, 0.0125, 0.025, 0.05, 0.1, 0.2, 0.4, 0.8, 1.0]}
    score_board = dict()
    index = 0
    min_fraction_to_keep = np.linspace(0, 0.5, 6)

    overall_ppv = []
    overall_recall = []
    overall_best_param = []

    # try parameters, pick the best combo
    for min_frac in tqdm(min_fraction_to_keep):
        min_frac_ppv = []
        min_frac_recall = []
        params = []
        for n in param_grid['num_neighbors']:
            for t in param_grid['t']:
                param_ppv = []
                param_recall = []
                no_fp_and_fn = 0
                for test_case in test_cases:
                    software_ibd_length, cohort_candidates, true_cohort_ibd_length, true_cohort_members = build_dataset_for_gridsearch(test_case)
                    knn_dict = build_knn_dict(cohort_candidates, single_edge_dict)

                    accepted_cohort_members = [x[0] for x in mean_dist_cohort(knn_dict, n, t, min_frac)]  
                    accepted_cohort_members = set(accepted_cohort_members)

                    true_cohort_members = set(true_cohort_members)

                    false_positives, true_positives, false_negatives = calculate_fp_tp_fn(accepted_cohort_members, true_cohort_members)

                    if (true_positives + false_positives) != 0:
                        ppv = true_positives/(true_positives + false_positives)
                        recall = true_positives/(true_positives + false_negatives)
                        param_ppv.append(ppv)
                        param_recall.append(recall)
                    else:
                        no_fp_and_fn +=1
                # calculate average!
                params.append((n, t, no_fp_and_fn))
                min_frac_ppv.append(np.mean(param_ppv))
                min_frac_recall.append(np.mean(param_recall))

        # calculate best 
        
        best_min_frac_idx = np.argmax(min_frac_ppv)
        overall_ppv.append(min_frac_ppv[best_min_frac_idx])
        overall_recall.append(min_frac_recall[best_min_frac_idx])
        overall_best_param.append(params[best_min_frac_idx])

    print(overall_ppv)
    logging.info(str(overall_ppv))
    print(overall_recall)
    logging.info(str(overall_recall))
    print(overall_best_param)
    logging.info(str(overall_best_param))
    print("-----")



