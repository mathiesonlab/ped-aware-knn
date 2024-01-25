from collections import deque
# class for one individual

class Individual:
    def __init__(self, indv_id, father_id, mother_id, children_ids = set()):
        self.indv_id = indv_id
        self.father_id = father_id
        self.mother_id = mother_id
        self.children_ids = children_ids

    def __hash__(self):
        return hash(self.indv_id)
    
    def __str__(self):
        return "{0} with father {1} and mother {2} and childrens {3}".format(
            self.indv_id, str(self.father_id), str(self.mother_id), 
            ", ".join(str(item) for item in self.children_ids))

# Pedigree operations

def read_ped_file(ped_filename):
    pedigree_dict = {}
    with open(ped_filename, 'r') as f:
        for row in f:
            info = row.strip().split()
            indv_id, father_id, mother_id = info[:3]

            if not indv_id in pedigree_dict:
                pedigree_dict[indv_id] = Individual(indv_id, father_id, mother_id)
            else:
                pedigree_dict[indv_id].mother_id = mother_id
                pedigree_dict[indv_id].father_id = father_id

            if not mother_id in pedigree_dict:
                pedigree_dict[mother_id] = Individual(mother_id, None, None, {indv_id})
            else:
                mother = pedigree_dict[mother_id]
                mother.children_ids.add(indv_id)
            
            if not father_id in pedigree_dict:
                pedigree_dict[father_id] = Individual(father_id, None, None, {indv_id})
            else:
                father = pedigree_dict[father_id]
                father.children_ids.add(indv_id)
            
    return pedigree_dict


def find_min_path(pedigree_dict, indv1, indv2):
    queue = deque()
    visited = set()
    ancestors = set()
    
    queue.append(indv1)
    ancestors.add(indv1)
    num_generation = 0

    while queue:
        current_generation = len(queue)
        while current_generation > 0:
            relative_id = queue.popleft()
            current_generation -=1
            if relative_id == indv2:
                return num_generation
            elif relative_id in visited:
                continue
            else:
                visited.add(relative_id)

                # add parents onto the pedigree queue
                r = pedigree_dict[relative_id]
                
                if relative_id in ancestors and relative_id != "0":
                    if r.father_id != "0":
                        queue.append(r.father_id)
                        ancestors.add(r.father_id)

                    if r.mother_id != "0":
                        queue.append(r.mother_id)
                        ancestors.add(r.mother_id)
                
                for child_id in r.children_ids:
                    queue.append(child_id)
        num_generation +=1
    return 0

