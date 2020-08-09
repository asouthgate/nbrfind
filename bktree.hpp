#include "uk2.hpp"
#include <vector>
#include <string>
#include <utility>

struct SearchResult {
    int d = -1;
    int NM = -1;
    int NN = -1;
    double p = -1;    
    int recind = -1;
    SearchResult(int d_, int NM_, int NN_, double p_, int recind_) (d=d_,NM=NM_,NN=NN_,p=p_,recind=recind_) {}
}

class Node {
    public:
        void add(const std::string& s);
        bool is_terminal();
        int recind = -1;
    private:
}

class BKTree {
    public:
        BKTree();
        void build(const std::string& ref_fasta_fname);
        std::list<Node*> query_add(const std::string& s, int k=5, float eps=0.95);
        std::vector<SearchResult> search(const std::string& s, int k=5, float eps=0.95);
    private:
        Node* root_ptr; 
        std::vector<std::pair<std::string, std::string>> refs;
}

