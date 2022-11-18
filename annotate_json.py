'''
This is adapted from https://github.com/jmcbroome/automate-lineages-prototype for application to 
Auspice formatted JSON with reduced features and minimal additional dependencies.
'''

import sys
import json
import argparse
from queue import SimpleQueue

def argparser():
    parser = argparse.ArgumentParser(description="Simple implementation of the genotype representation metric for automated lineage designation for arbitrary Nextstrain JSON.")
    parser.add_argument("-i","--input",help="Name of an input JSON.",required=True)
    parser.add_argument("-o","--output",help="Name of an output annotated JSON.",required=True)
    parser.add_argument("-f","--floor",type=int,default=0,help="Set a minimum total value to annotate a lineage.")
    parser.add_argument("-s","--size",type=int,default=0,help="Set a minimum number of samples to annotate a lineage.")
    parser.add_argument("-d","--distinction",type=int,default=0,help="Set a minimum number of mutations separating a new lineage label with its parent.")
    parser.add_argument("-c","--cutoff",type=float,default=1,help="Proportion of samples that must be labeled on each level.")
    parser.add_argument("-m","--missense",action='store_true',default=False,help="Use to only consider amino-acid altering mutations.")
    parser.add_argument("-g","--gene",default=None,help="Only consider missense mutations within a specific gene. Sets -m")
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
        if c.id not in banned:
            cscore = evaluate_candidate(anid, c.id, sum_and_count, dist_to_root, minimum_size, minimum_distinction)
            # if cscore > 50:
                # print(c.id, cscore)
            if cscore > 0:
                good_candidates.append((cscore,c))
    # if len(good_candidates) > 0:
        # print("CANDIDATEINFO", len(good_candidates))
    if len(good_candidates) == 0:
        return (0,None)
    return max(good_candidates, key=lambda x: x[0])

def update_json(ijd, labels, levels=1):
    for l in range(1,levels):
        ijd['meta']['colorings'].append({"key":"autolin_level_"+str(l),"title":"autolin_level_"+str(l),"type":"categorical"})
    treed = ijd['tree']
    def traverse(cnd):
        if "name" in cnd.keys():
            flabel = labels.get(cnd['name'],'L')
            for l in range(1,levels):
                stripped = ".".join(flabel.split(".")[:l+1])
                cnd['node_attrs']['autolin_level_'+str(l)] = {'value':stripped}
        for nd in cnd.get("children",[]):
            traverse(nd)
    traverse(treed)
    ijd['treed'] = treed
    return ijd

class TreeNode:
    def __init__(self, nid, parent=None, mutations=[]):
        self.id = nid
        self.mutations = mutations
        self.children = []
        self.parent = parent

    def add_child(self, obj):
        self.children.append(obj)

    def add_mutation(self, obj):
        self.mutations.append(obj)

    def is_leaf(self):
        return (len(self.children) == 0)

    def __str__(self, level=0):
        ret = "\t"*level+repr(self.id)+"\n"
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
        
    def __loader(self, jd, cnode, aa=False, gene=None):
        muinfo = jd['branch_attrs']['mutations']
        global id_counter
        if not aa and gene == None:
            if 'nuc' in muinfo.keys():
                for m in muinfo['nuc']:
                    cnode.add_mutation(m)
        else:
            for g, aav in muinfo.items():
                if g != 'nuc':
                    if gene == None or g == gene:
                        for aa in aav:
                            cnode.add_mutation(aa)
        for child in jd.get("children",[]):
            if 'name' in child.keys() and 'children' not in child.keys():
                new_nid = child['name']
            else:
                new_nid = 'node_' + str(id_counter)
            id_counter += 1
            child_node = self.__loader(child, TreeNode(new_nid, parent=cnode), aa, gene)
            cnode.add_child(child_node)
            self.nodes[child_node.id] = child_node
        return cnode

    def load_from_dict(self, jd, nid_ccount = 1, aa = False, gene = None):
        global id_counter
        id_counter = nid_ccount
        cnode = self.root
        self.__loader(jd, cnode, aa, gene)
        return self

    def get_node(self, nid):
        return self.nodes.get(nid, None)

    def parsimony_score(self):
        return sum([len(n.mutations) for n in self.nodes.values()])

    def rsearch(self, node):
        cp = node
        path = [cp.id]
        while cp.parent != None:
            path.append(cp.parent.id)
            cp = cp.parent
        return path

    def breadth_first_expansion(self, cnode=None, reverse=False):
        if cnode == None:
            cnode = self.root
        bfs = []
        remaining = SimpleQueue()
        remaining.put(cnode)
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
        return leaf_ids

    def __str__(self):
        return self.root.__str__()

def pipeline(ijson, ojson, floor=0, size=0, distinction=0, cutoff=1, missense=False, gene=None):
    with open(ijson) as inf:
        ijd = json.load(inf)
    t = Tree().load_from_dict(ijd['tree'], 1, missense, gene)
    print(f"Loaded tree successfully; parsimony score {t.parsimony_score()}.",file=sys.stderr)
    annotes = {'L':t.root.id}
    outer_annotes = annotes
    level = 1
    all_labels = {}
    while True:
        new_annotes = {}
        used_nodes = set()
        for ann,nid in outer_annotes.items():
            serial = 0
            rbfs = t.breadth_first_expansion(t.get_node(nid), True)
            # print(f"Breadth first expansion complete, size {len(rbfs)}.",file=sys.stderr)
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
                for anc in t.rsearch(best_node):
                    used_nodes.add(anc)
                # print(used_nodes)
                new_annotes[newname] = best_node.id
                leaves = t.get_leaves_ids(best_node)
                # print(f"Leaf fetching complete- {len(leaves)} leaves.")
                for l in leaves:
                    labeled.add(l)
                    #overrwite an existing higher-level label if it exists
                    #beacuse each lineage label name contains its ancestral lineage labels as well.
                    all_labels[l] = newname
                # print("LABELINFO", len(labeled), leaf_count, cutoff)
                if len(labeled) >= leaf_count * cutoff:
                    break
                serial += 1
                print(f"Annotation {newname} generated for node {best_node.id}.")
        if len(new_annotes) == 0:
            break
        else:
            annotes.update(new_annotes)
            outer_annotes = new_annotes
            level += 1

    print(f"Total samples labeled: {len(all_labels)}\nTotal labels generated: {len(annotes)}")
    njd = update_json(ijd, all_labels, level)
    with open(ojson,'w+') as of:
        json.dump(njd,of)
    
def main():
    args = argparser()
    pipeline(args.input,args.output,args.floor,args.size,args.distinction,args.cutoff,args.missense,args.gene)

if __name__ == "__main__":
    main()