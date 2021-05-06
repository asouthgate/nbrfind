#include <vector>
#include <string>
#include <utility>
#include <iostream>
#include "sparseMEM_src/sparseSA.hpp"
#include "sdata.hpp"
#include "cluster.hpp"

const size_t MIN_LENGTH = 29000;
const size_t MAX_LENGTH = 31000;
const size_t MAX_ROW_SIZE = 2* MAX_LENGTH + 1;

struct unalignedSegment {
    int i_start = -1;
    int i_end = -1;
    int j_start = -1;
    int j_end = -1; 
};

struct MoveResult {
    MoveResult(int maxi_, int prev_NM_, int prev_NN_) : maxi(maxi_), prev_NM(prev_NM_), prev_NN(prev_NN_) {}
    MoveResult() {}
    // Index of the optimal move
    int maxi = -1;
    // Defining the previous number of N-X mismatches and number of N-N matches
    int prev_NM = -1;
    int prev_NN = -1;
    // Success - unnecessary
//    bool success = true;
};

std::vector<unalignedSegment> extract_unmatched_segments(std::vector<match_t> matches, std::string S, int m, int n);


class DistCalculator {
    public:
        DistCalculator() {}
        bool calculate_dist_sd(std::string s1, std::string s2, StateData& sd, int snpmax, int slide_threshold, bool freeze=false);
        void query_samples_against_refs(std::string sample_fasta_fname, std::string ref_fasta_fname, int k=5, int max_slide=100, double epsilon = 0.5);
        std::vector<Cluster> get_clusters(const std::vector<std::pair<std::string,std::string>>& refs, int k, double epsilon);

        std::tuple<int,int,int,int> align2seqs_mem(StateData& sd, sparseSA& spsa, std::string& p1, std::string& p2, int k, int max_slide);

        std::tuple<int,int,int,int> align2seqs(StateData& sd, std::string& p1, std::string& p2, int k, int max_slide);

    private:
};
