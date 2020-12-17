#include<iostream>
#include<vector>
#include<string>
#include<assert.h>
#include"uk2.hpp"
#include"seqio.hpp"
#include "sparseMEM_src/sparseSA.hpp"


int main(int argc, char* argv[]) {
    std::cerr << "testing" << std::endl; 
    std::string S = "NATATTNNATANANAAAATTTTTN";
    std::pair<std::vector<int>, std::string> res = strip_Ns(S);
    for (int j = 0; j < res.first.size(); ++j) {
        std::cerr << j << " " << res.second[j] << " " << res.first[j] << std::endl;
        if (res.first[j] >= 0) assert(res.second[j] == S[res.first[j]]);
    }

    S = "NATATTNNATANANAAAATTTTTNN";
    res = strip_Ns(S);
    for (int j = 0; j < res.first.size(); ++j) {
        std::cerr << j << " " << res.second[j] << " " << res.first[j] << std::endl;
        if (res.first[j] >= 0) assert(res.second[j] == S[res.first[j]]);
    }

    S = "ATATTNNATANANAAAATTTTTNN";
    res = strip_Ns(S);
    for (int j = 0; j < res.first.size(); ++j) {
        std::cerr << j << " " << res.second[j] << " " << res.first[j] << std::endl;
        if (res.first[j] >= 0) assert(res.second[j] == S[res.first[j]]);
    }

    S = "ATATTNNATANANAAAATTTTT";
    res = strip_Ns(S);
    for (int j = 0; j < res.first.size(); ++j) {
        std::cerr << j << " " << res.second[j] << " " << res.first[j] << std::endl;
        if (res.first[j] >= 0) assert(res.second[j] == S[res.first[j]]);
    }


    std::vector<std::pair<std::string, std::string>> recs1 = read_fasta(argv[1], 29000,32000);
    std::vector<std::pair<std::string, std::string>> recs2 = read_fasta(argv[2], 29000,32000);
    std::string rec1 = recs1[0].second;

    // Sparse SA 
    std::vector<long> startposs = {0};
    std::vector<std::string> descriptions = {"foo"};
    sparseSA spsa = sparseSA(rec1, descriptions, startposs, false, 1);
    for (int ri = 0; ri < recs2.size(); ++ri) {
        std::string rec2 = recs2[ri].second;
        std::vector<match_t> matches;
        spsa.MUM(rec2, matches, 20, false);
        std::vector<unalignedSegment> uas = extract_unmatched_segments(matches, rec1.length(), rec2.length());

    }

    return 0;
}
