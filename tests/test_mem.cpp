#include<iostream>
#include<vector>
#include<string>
#include<assert.h>
#include<tuple>
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

    S = "NNNNNNATATTNNATANANAAAATTNNNNNNNTTT";
    res = strip_Ns(S);
    for (int j = 0; j < res.first.size(); ++j) {
        std::cerr << j << " " << res.second[j] << " " << res.first[j] << std::endl;
        if (res.first[j] >= 0) assert(res.second[j] == S[res.first[j]]);
    }

    std::cerr << "Testing short strings..." << std::endl;

    std::vector<std::string> short_strings = {"ATATTNAAATTNANANAATTTGTTN", "GTCTTTTNAAATTNANANAA", "ATATTNAAATTNANGTTTTTTAGAATTTTTTTTTTTNNNNNNNNNATTTGNNNTTNNNNNNNNN"};
    DistCalculator dc = DistCalculator();

    StateData short_sd(200 + 1);
    std::vector<long> startposs = {0};
    std::vector<std::string> descriptions = {"foo"};
    sparseSA short_spsa = sparseSA(S, descriptions, startposs, false, 1);

    for (auto& rec2 : short_strings) {
        std::cerr << "running short comparison..." << std::endl;
        std::tuple<int,int,int,int> seeded_alignment = dc.align2seqs_mem(short_sd, short_spsa, S, rec2, 1000, 50000);
        std::tuple<int,int,int,int> alignment = dc.align2seqs(short_sd, S, rec2, 1000, 50000);
        std::cerr << get<0>(seeded_alignment) << " " << get<0>(alignment) << std::endl;
        assert(get<0>(seeded_alignment) == get<0>(alignment));
    }

    std::cerr << "short strings passed!" << std::endl;

    std::cerr << "reading in seqs.." << std::endl;
    std::vector<std::pair<std::string, std::string>> recs1 = read_fasta(argv[1], 29000,32000);
    std::vector<std::pair<std::string, std::string>> recs2 = read_fasta(argv[2], 29000,32000);
    std::string rec1 = recs1[0].second;

    size_t max_seq_len = 0;
    for (int si = 0; si < recs2.size(); ++si) {
        max_seq_len = std::max(max_seq_len, recs2[si].second.length());
    }
    for (int si = 0; si < recs1.size(); ++si) {
        max_seq_len = std::max(max_seq_len, recs1[si].second.length());
    }

    std::cerr << "testing unaligned segments.." << std::endl;
    // Sparse SA 
    sparseSA spsa = sparseSA(rec1, descriptions, startposs, false, 1);
    for (int ri = 0; ri < recs2.size(); ++ri) {
        std::string rec2 = recs2[ri].second;
        std::vector<match_t> matches;
        spsa.MUM(rec2, matches, 100, false);
        std::vector<unalignedSegment> uas = extract_unmatched_segments(matches, rec1, rec1.length(), rec2.length());
        int i_prev_end = 0;
        int j_prev_end = 0;
        for (auto& u : uas) {
            std::cerr << u.i_start << " " << u.j_start << std::endl;
            assert(u.i_start >= i_prev_end);
            assert(u.j_start >= j_prev_end);
            j_prev_end = u.j_end;
            i_prev_end = u.i_end;
        }
    }
    std::cerr << "Test passed, no mums overlap" << std::endl;

    std::cerr << "final test; real sequence comparison" << std::endl;
    StateData sd(max_seq_len*2 + 1);
    for (auto& p : recs2) {
        std::string rec2 = p.second;
        std::cerr << "running comparison..." << std::endl;
        std::tuple<int,int,int,int> seeded_alignment = dc.align2seqs_mem(sd, spsa, rec1, rec2, 1000, 50000);
        std::tuple<int,int,int,int> alignment = dc.align2seqs(sd, rec1, rec2, 1000, 50000);
        std::cerr << "should be equal... " << get<0>(seeded_alignment) << " " << get<0>(alignment) << std::endl;
        assert(get<0>(seeded_alignment) == get<0>(alignment));
    }

    return 0;
}
