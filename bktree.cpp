#include "uk2.hpp"
#include <vector>
#include <string>

SearchResult cal_p(std::string x, std::string y) {
    int rowsize = x.length() + y.length() + 1;
    init_state_array(state_arr, rowsize);
    init_state_triple(state_triple, x.length(), y.length());
    int mm_ind = 5*rowsize + y.length();
    int NM_ind = 8*rowsize + y.length();
    int NN_ind = 11*rowsize + y.length();
    calculate_dist(x, y, state_triple, state_arr, rowsize, 5, 100, false);
    int d = state_arr[mm_ind];
    int NM = state_arr[NM_ind];
    int NN = state_arr[NN_ind];
    double p = betap(d, NM, NN, k);
    return SearchResult(d,NM,NN,p);
}

bool Node::is_terminal() {
    if (children.size() > 0) { return true;}
    return false;
}

std::list<Node*> BKTree::query_add(const std::string& s, int k=5, float eps=0.95) {
    std::list<Node*> node_stack;
    // add all children of the root
    for (auto & c: root_ptr->children) {
        node_stack.push_back(c);
    }
    std::list<Node*> results;
    while (node_stack.size() > 0) {
        Node* currn_ptr = node_stack.pop_back();
        // Query
        int matching_children = 0;
        for (auto & c : currn_ptr->children) {
            if (cal_p(s, refs[c->recind].second)) {
                if (c->terminal) {
                    results.push_back(c);
                    matching_children += 1;
                }
                else {
                    node_stack.push_back(c);
                }
            }
        }
        if (matchind_children == 0) results.push_back(currn_ptr);
        

    }
    return results;
}

void BKTree:build(const std::string& ref_fasta_fname, int k, float eps) {
    refs = read_fasta(ref_fasta_fname);
    std::string x0 = refs[0].second;
    for (auto& p2 : refs) {
        std::list<Node*> res = query_add(p2);
        for (Node* n : res) {
            n.add(p2)
        }
    }
}

