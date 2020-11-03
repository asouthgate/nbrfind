#include <vector>
#include <string>
#include <utility>
#include <iostream>
#include "sdata.hpp"
#include "cluster.hpp"

const size_t MIN_LENGTH = 29000;
const size_t MAX_LENGTH = 31000;
const size_t MAX_ROW_SIZE = 2* MAX_LENGTH + 1;

struct MoveResult {
    MoveResult(int maxi_, int prev_NM_, int prev_NN_) : maxi(maxi_), prev_NM(prev_NM_), prev_NN(prev_NN_) {}
    MoveResult() {}
    int maxi = -1;
    int prev_NM = -1;
    int prev_NN = -1;
    bool success = true;
};

class DistCalculator {

    public:
        DistCalculator() {}
        bool calculate_dist_sd(std::string s1, std::string s2, StateData& sd, int snpmax, int slide_threshold, bool freeze=false);
        void query_samples_against_refs(std::string sample_fasta_fname, std::string ref_fasta_fname, int k=5, double epsilon = 0.5);
        std::vector<Cluster> get_clusters(const std::vector<std::pair<std::string,std::string>>& refs, int k, double epsilon);

    private:
};
