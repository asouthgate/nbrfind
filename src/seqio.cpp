#include <string>
#include <iostream>
#include <cstring>
#include <fstream>
#include <algorithm>
#include <vector>
#include "seqio.hpp"

std::vector<std::pair<std::string, std::string>> read_fasta(std::string fasta_fname, int min_len, int max_len) {
    std::cerr << "Reading " << fasta_fname << std::endl;
    std::vector<std::pair<std::string,std::string>> records;
    std::string curr_str = "";
    std::string curr_h = "";
    std::string l;
    int rec_counter = 0;
    int c = 0;
    std::ifstream fasta_file(fasta_fname);
    while (std::getline(fasta_file,l)) {
        l.erase(std::remove(l.begin(), l.end(), '\n'), l.end());
        l.erase(std::remove(l.begin(), l.end(), '\r'), l.end());
        if (l[0] == '>') {
            if (c > 0) {
//                std::transform(curr_str.begin(), curr_str.end(), curr_str.begin(),
//               [](unsigned char c){ return std::toupper(c); });
//                convert_nonstandard_to_N(curr_str);
//                bool cds = trim_to_cds(curr_str);
//                if (cds && filter(curr_str, min_len, max_len)) {
                records.push_back(std::pair<std::string,std::string>(curr_h, curr_str));
//                }
                rec_counter += 1;
            }
            curr_h = l;
            curr_str = "";
        }
        else {
            curr_str = curr_str + l;
        }
        c += 1;
    }
    rec_counter += 1;
//    std::transform(curr_str.begin(), curr_str.end(), curr_str.begin(),
//           [](unsigned char c){ return std::toupper(c); });
//    convert_nonstandard_to_N(curr_str);
//    bool cds = trim_to_cds(curr_str);
//    if (cds && filter(curr_str, min_len, max_len)) {
    std::pair<std::string,std::string> record(curr_h, curr_str);
    records.push_back(record);
//    }
    std::cerr << "Excluded " << rec_counter-records.size() << " of " << records.size() << " records due to filter..." << std::endl;
    return records;
}



bool trim_to_cds(std::string& s) {
    std::size_t found5 = s.find(five_prime);
    if (found5!=std::string::npos) {
        s = s.substr(found5);
    }
    else { std::cerr << "failed to find 3'" << std::endl; return false; }
    std::size_t found3 = s.find(three_prime);
    if (found3!=std::string::npos) {
        s = s.substr(0,found3+three_prime.length());
    }
    else { std::cerr << "failed to find 3'" << std::endl; return false; }
    return true;
}

void convert_nonstandard_to_N(std::string & s) {
    for (int j = 0; j < s.size(); ++j) {
        if (s[j] == 'A') {}
        else if (s[j] == 'C') {}
        else if (s[j] == 'G') {}
        else if (s[j] == 'T') {}
        else s[j] = 'N';
    }
}

bool filter(const std::string& s, int min_len, int max_len) {
    if (s.length() < min_len) return false;
    if (s.length() > max_len) return false;
    size_t n = std::count(s.begin(), s.end(), 'N');
    float perc_gaps = n/ (float) s.length();
    if (perc_gaps > 0.2) {
        std::cerr << "too many gaps'" << std::endl; return false;
    }
    return true;
}


