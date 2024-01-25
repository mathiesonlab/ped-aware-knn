from pedigree import find_min_path
from statistics import mean

def build_knn_dict(cohort, pedigree_dict):
    knn_dict = dict()
    for person in cohort:
        if knn_dict.get(person) is None:
            knn_dict[person] = dict()
        for another_person in cohort:
            if another_person == person or knn_dict[person].get(another_person) is not None:
                pass
            else:
                knn_dict[person][another_person] = find_min_path(pedigree_dict, person, another_person)
    return knn_dict

def mean_dist_cohort(knn_dict, num_neighbors, t, min_fraction_to_keep):
    indiv_dict_with_mean_knn_num_edge = dict()
    for key, value in knn_dict.items():
        sorted_min_edge = sorted(list(value.values()))
        if len(sorted_min_edge) > 0:
            mean_num_edge = mean(sorted_min_edge[:num_neighbors])
            indiv_dict_with_mean_knn_num_edge[key] = mean_num_edge
    ordered_list_of_indivs_with_mean_knn_num_edge = [(k, v) for k, v in sorted(indiv_dict_with_mean_knn_num_edge.items(), key=lambda item: item[1])]
    largest_increase_in_mean_knn_num_edge = 0
    list_of_increase_mean_knn_num_edge_for_adj_indiv = []
    for i in range(1, len(ordered_list_of_indivs_with_mean_knn_num_edge)):
        increase_in_mean_knn_num_edge = ordered_list_of_indivs_with_mean_knn_num_edge[i][1] - ordered_list_of_indivs_with_mean_knn_num_edge[i-1][1]
        list_of_increase_mean_knn_num_edge_for_adj_indiv.append(increase_in_mean_knn_num_edge)
        if increase_in_mean_knn_num_edge > largest_increase_in_mean_knn_num_edge:
            largest_increase_in_mean_knn_num_edge = increase_in_mean_knn_num_edge
    threshold = largest_increase_in_mean_knn_num_edge * t
    critical_index = 0
    for i in range(len(list_of_increase_mean_knn_num_edge_for_adj_indiv)):
        if list_of_increase_mean_knn_num_edge_for_adj_indiv[i] >= threshold:
            critical_index = i
            break
    
    final_index = max(round(min_fraction_to_keep * len(ordered_list_of_indivs_with_mean_knn_num_edge)), 1, critical_index)
    return ordered_list_of_indivs_with_mean_knn_num_edge[:final_index+1] # check +- 1 for critical index here



