#include<vector>
#include<string>
#include<utility>

const size_t MIN_LENGTH = 29000;
const size_t MAX_LENGTH = 31000;
const size_t MAX_ROW_SIZE = 2* MAX_LENGTH + 1;

class DistCalculator {
    public:
        DistCalculator() {}
        bool calculate_dist(std::string s1, std::string s2, int* state_triple, int* state_arr, int rowsize, int snpmax, int slide_threshold, bool freeze=false);
        void query_samples_against_refs(std::string sample_fasta_fname, std::string ref_fasta_fname);
    private:
        std::vector<std::pair<std::string, std::string>> read_fasta(std::string fasta_fname);
        void init_state_array(int* state_arr, int rowsize);
        void init_state_triple(int* state_triple, int len1, int len2);
        double betap(int d, int M, int N, int k, int a=1, int b=1);
};
