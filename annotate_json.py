import sys
import json
import argparse
from queue import SimpleQueue

def argparser():
    parser = argparse.ArgumentParser(description="Simple implementation of the genotype representation metric for automated lineage designation for arbitrary Nextstrain JSON.")
    parser.add_argument("-i","--input",help="Name of an input JSON.",required=True)
    parser.add_argument("-o","--output",help="Name of an output annotated JSON.",required=True)
    parser.add_argument("-m","--missense",help="Consider only amino acid substituting mutations.")
    parser.add_argument("-f","--floor",help="Set a minimum total value to annotate a lineage.")
    parser.add_argument("-s","--size",help="Set a minimum number of samples to annotate a lineage.")
    parser.add_argument("-d","--distinction",help="Set a minimum number of mutations separating a new lineage label with its parent.")
    return parser.parse_args() 

def dists_to_root(node):
    #nodes must be a dict that gets updated on each recursion
    #gives back a dict with all nodes and their respective dist from root
    #initalize this with our starting node at 0, its our "root" whether its the actual tree root or not
    nodes = {node.id:0}
    def recursive_dists_to_roots(snode):
        bweight = nodes[snode.id]
        for child in snode.children:
            dist = bweight + len(child.mutations) 
            nodes[child.id] = dist
            recursive_dists_to_roots(child)
    recursive_dists_to_roots(node)
    return nodes

def get_sum_and_count(rbfs, ignore = set()):
    # node sum stored in first index and node count stored in second index of each dict entry
    sum_and_count_dict = {}
    leaf_count = 0
    for node in rbfs:
        if node.is_leaf():
            leaf_count += 1
            if node.id not in ignore:
                sum_and_count_dict[node.id] = (len(node.mutations), 1)
        else:
            total_count = 0
            total_sum = 0
            for child in node.children:
                sumtc = sum_and_count_dict.get(child.id, None)
                if sumtc == None:
                    continue
                total_count += sumtc[1]
                total_sum += sumtc[0]
            if total_count > 0:
                #total path length is computed as the total path lengths to each child plus the length of the current node TIMES the number of samples.
                #this is because total path length is not the same as tree parsimony- some mutations are part of many sample paths
                #for a given sample to its parent, the total path length is just the number of mutations (as computed above)
                #but for an internal node with two leaf children's path length with respect to its parent, 
                #its equal to the sum of the two child's path lengths plus 2 times its mutations, since those mutations are shared among 2 samples
                #this logic applies as we move further up the tree.
                sum_and_count_dict[node.id] = (total_sum + len(node.mutations) * total_count, total_count)
    return sum_and_count_dict, leaf_count

def evaluate_candidate(a, nid, sum_and_counts, dist_to_root, minimum_size=0,minimum_distinction=0):
    """Evaluate a candidate branch as a putative sublineage.

    Args:
        t (MATree): The tree.   
        a (str): The parent lineage annotation node.
        nid (str): The node id of the candidate branch.
    """
    node_sum, node_count = sum_and_counts.get(nid,[0,0])
    if node_count <= minimum_size:
        return 0
    if node_sum == 0 or node_count <= 0:
        return 0
    candidate_to_parent = dist_to_root[nid] - dist_to_root[a]
    if candidate_to_parent < minimum_distinction:
        return 0
    mean_distances = node_sum/node_count
    if (mean_distances + candidate_to_parent) == 0:   #avoid divide by 0
        candidate_value = 0
    else:
        # print("DEBUG: {} {} {} {}".format(node_count, max([(node_count-minimum_size+1),0]), candidate_to_parent, max([candidate_to_parent-minimum_distinction+1,0])))
        # candidate_value = max([(node_count-minimum_size+1),0]) * (max([candidate_to_parent-minimum_distinction+1,0])) / (mean_distances + candidate_to_parent)
        candidate_value = node_count * candidate_to_parent / (mean_distances + candidate_to_parent)
    return candidate_value

def evaluate_lineage(t, dist_to_root, anid, candidates, sum_and_count, minimum_size = 0, minimum_distinction = 0, banned = set()):
    good_candidates = []
    for c in candidates:
        if not c.is_leaf() and c.id not in banned:
            cscore = evaluate_candidate(anid, c.id, sum_and_count, dist_to_root, minimum_size, minimum_distinction)
            if cscore > 0:
                good_candidates.append((cscore,c))
    if len(good_candidates) == 0:
        return (0,None)
    return max(good_candidates, key=lambda x: x[0])

def update_json(ijd, labels):
    treed = ijd['tree']
    def traverse(cnd):
        if "name" in cnd.keys():
            cnd['node_attrs']['autolin'] = {'value':labels.get(cnd['name'],'None')}
        for nd in treed['children']:
            traverse(nd)
    traverse(treed)
    return treed

class TreeNode:
    def __init__(self, nid, mutations=[]):
        self.id = nid
        self.mutations = mutations
        self.children = []

    def add_child(self, obj):
        self.children.append(obj)

    def add_mutation(self, obj):
        self.mutations.append(obj)

    def is_leaf(self):
        return (len(self.children) == 0)

    def __str__(self, level=0):
        ret = "\t"*level+repr(self.nid)+"\n"
        for child in self.children:
            ret += child.__str__(level+1)
        return ret

    def __repr__(self):
        return "\n".join(["id: "+self.id,"# of mutations: "+str(len(self.mutations)),"# of children: "+str(len(self.children))])

class Tree:
    '''
    Minimalist tree class supporting the application of the genotype representation heuristic for lineage nomenclature creation.
    '''
    def __init__(self):
        self.root = TreeNode('node_0')
        self.nodes = {'node_0':self.root}
        
    def __loader(self, jd, cnode):
        for m in jd['branch_attrs']['mutations']['nuc']:
            cnode.add_mutation(m)
        for child in jd['tree']['children']:
            if 'name' in child.keys():
                new_nid = child['name']
            else:
                new_nid = 'node_' + str(id_counter)
            id_counter += 1
            child_node = self.__loader(child, TreeNode(new_nid))
            cnode.add_child(child_node)
            self.nodes[child_node.id] = child_node
        return cnode

    def load_from_dict(self, jd, nid_ccount = 1):
        global id_counter
        id_counter = nid_ccount
        cnode = self.root
        self.__loader(jd, cnode)

    def get_node(self, nid):
        return self.nodes.get(nid, None)

    def breadth_first_expansion(self, cnode=None, reverse=False):
        if cnode == None:
            cnode = self.root
        bfs = []
        remaining = SimpleQueue()
        while not remaining.empty():
            node = remaining.get()
            bfs.append(node)
            for c in node.children:
                remaining.put(c)
        if reverse:
            bfs.reverse()
        return bfs
        
    def get_leaves_ids(self,cnode=None):
        allnodes = self.breadth_first_expansion(cnode,reverse=True)
        leaf_ids = []
        for n in allnodes:
            if n.is_leaf():
                leaf_ids.append(n.id)
            else:
                #in reverse breadth-first, all leaves should appear first. 
                #therefore, when leaves run out, we are done.
                break
        return leaf_ids

    def __str__(self):
        return self.root.__str__()

def pipeline(ijson, ojson, floor=0, size=0, distinction=0, cutoff=1):
    ijd = json.load(ijson)
    t = Tree().load_from_dict(ijd['tree'])
    annotes = {'L':'node_0'}
    outer_annotes = annotes
    level = 1
    all_labels = {}
    while True:
        new_annotes = {}
        used_nodes = set()
        for ann,nid in outer_annotes.items():
            serial = 0
            rbfs = t.breadth_first_expansion(nid, True)
            parent_leaf_count = len([n for n in rbfs if n.is_leaf()])
            if parent_leaf_count == 0:
                continue
            labeled = set()
            dist_root = dists_to_root(t.get_node(nid)) #needs the node object, not just the name
            while True:
                scdict, leaf_count = get_sum_and_count(rbfs, ignore = labeled)
                best_score, best_node = evaluate_lineage(t, dist_root, nid, rbfs, scdict, size, distinction, used_nodes)
                if best_score <= floor:
                    break
                newname = ann + "." + str(serial)
                for anc in t.rsearch(best_node.id,True):
                    used_nodes.add(anc.id)
                new_annotes[newname] = best_node.id
                leaves = t.get_leaves_ids(best_node.id)
                for l in leaves:
                    labeled.add(l)
                    #overrwite an existing higher-level label if it exists
                    #beacuse each lineage label name contains its ancestral lineage labels as well.
                    all_labels[l] = newname

                if len(labeled) >= leaf_count * cutoff:
                    break
                serial += 1
        if len(new_annotes) == 0:
            break
        else:
            annotes.update(new_annotes)
            outer_annotes = new_annotes
            level += 1

    njd = update_json(ijd, all_labels)
    with open(ojson,'w+') as of:
        json.dump(njd,of)
    
def main():
    args = argparser()
    pipeline(args.input,args.output,args.missense,args.floor,args.size,args.distinction)

if __name__ == "__main__":
    main()