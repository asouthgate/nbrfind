#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1
#include <cmath>
#include <assert.h>
#include <math.h>
#include <algorithm>
#include <string>
#include <iostream>
#include <cstring>
#include <fstream>
#include "uk2.hpp"

int slide(int d, int i, const std::string& s1, const std::string& s2) {
    // Fastest way to optimize this?
    for (int q = i; q < std::min(s1.length(), s2.length()-d); ++q) {
        if ((s1[q] != s2[q+d]) && (s1[q] != 'N') && (s2[q+d] != 'N')){ 
            return q;
        }
    }
    return s1.length();
}

void cerrarr(int* v, int rowsize) {
    for (int pj = 0; pj < rowsize; ++pj) {
        std::cerr << v[pj] << " ";
    }
    std::cerr << std::endl;
}

unsigned long long binom(unsigned int n, unsigned int k) {
    if (n == k || k == 0) return 1;
    return C(n-1, k-1) * n/k;
}

double betap(d,M,N,k,a=1,b=1) {
    if (d >= k) return 0;
    double sum = 0;
    for (int z = 0; z < k; ++z) {
        unsigned long long b = binom(N,z);
        long double numo = std::betal(z-a+d,N-z+b+M-d);
        long double deno = std::betal(a+d, b+M-d);
    }
} 

bool DistCalculator::calculate_dist(std::string s1, std::string s2, int* state_triple, int* state_arr, int rowsize, int snpmax, int slide_threshold, bool freeze) {
    int m = s1.length();
    int n = s2.length();
    // make sure frozen not true
    // init
    int h_start = state_triple[0];
    int lower_bound = state_triple[1];
    int upper_bound = state_triple[2];
    auto L0 = state_arr;
    auto L1 = state_arr + rowsize;
    auto L2 = state_arr + 2*rowsize;
    auto M0 = state_arr + 3*rowsize;
    auto M1 = state_arr + 4*rowsize;
    auto M2 = state_arr + 5*rowsize;
//    std::cerr << "RESUMING:" << "h=" << h_start << ",lb=" << lower_bound << ",ub=" << upper_bound << std::endl; 
//    std::cerr << "m n" << m << " " << n << std::endl;
//    std::cerr << "FREEZE?" << freeze << std::endl;
//    cerrarr(L0, rowsize);
//    cerrarr(L1, rowsize);
//    cerrarr(L2, rowsize);
//    std::cerr << std::endl;
//    cerrarr(M0, rowsize);
//    cerrarr(M1, rowsize);
//    cerrarr(M2, rowsize);

    if (freeze && n == 0) { return false; }

    int h = h_start;
    while (h < 2*(m+n)+1) {
        state_triple[0] = h;
//    for (int h = h_start; h < 2*(m+n)+1; ++h) {
        // printarr
//        std::cout <<  "h=" << h << " " << std::endl;
//        std::cout << std::max(-h/2, lower_bound) << " " << h << " " << std::min(h/2 + 1, upper_bound+1) << std::endl;
//        cerrarr(L0, rowsize);
//        cerrarr(L1, rowsize);
//        cerrarr(L2, rowsize);
//        std::cerr << std::endl;
//        cerrarr(M0, rowsize);
//        cerrarr(M1, rowsize);
//        cerrarr(M2, rowsize);
//        std::cerr << std::endl;
        int prev_lower_bound = lower_bound;
        int prev_upper_bound = upper_bound;
        for (int d = std::max(-(h/2), lower_bound); d < std::min((h/2)+1,upper_bound+1); ++d) {
 //           std::cout << "\t" << lower_bound << " " << d << " " << upper_bound << std::endl;
            int ld = m+d;
            int maxi = -1;
            int lmove = L0[ld-1];
            int matchmove = L1[ld] + 1;
            int prev_mm = M1[ld];
            int rmove = -1;
            if ((d > -h) && (ld-1 > 0)) {
                maxi = std::max(maxi, lmove);
                M2[ld] = M0[ld-1];
            }
            if ((d < h) && (ld+1 < m+n+1)) {
                rmove = L0[ld+1] + 1;
                maxi = std::max(maxi, rmove);
                if (rmove >= lmove) { M2[ld] = M0[ld+1]; }
            }
            if (matchmove >= maxi) {
                maxi = matchmove;
                if (maxi > 0) {
                    M2[ld] = prev_mm + 1;
                }
            }
//            std::cout << "maxi" << maxi << std::endl;
            L2[ld] = slide(d, maxi, s1, s2);
            int imax = m - std::max(0, d-(n-m));
            if (L2[ld] > imax) { L2[ld] = imax; }
            if (freeze && d >= n-m) {
                if (L2[ld] == imax) {
//                    std::cout << "freezing" << std::endl;
                    return true;
                }
            }
            assert (M2[ld] >= 0);
            if (d < 0 && L2[ld] >= m && d > lower_bound) {
                lower_bound = d; state_triple[1] = d;
            }
            if (d > 0 && L2[ld] >= n-d && d < upper_bound) {
                upper_bound = d-1; state_triple[2] = d-1;
            } 
        }
        // TODO: figure out fastest way to copy/clear the three arrays
//       cerrarr(L2);
//        std::cout << std::endl;
///        std::cout << "checking end" << L2[n] << "=? " << m << std::endl;
//        std::cout << "n:" << n << std::endl;
//        cerrarr(L0, rowsize);
//       cerrarr(L1, rowsize);
//        cerrarr(L2, rowsize);
//        std::cout << std::endl;
//        cerrarr(M0, rowsize);
//        cerrarr(M1, rowsize);
//        cerrarr(M2, rowsize);        
        if (L2[n] == m) {
//            std::cout << "ending early at h" << " " << h << std::endl;
            return true;
        }
//        L0 = L1;
//        std::cout << "memcopying" << std::endl;
        std::memmove(L0, L1, rowsize * sizeof(L0[0]));
        std::memmove(L1, L2, rowsize * sizeof(L0[0]));
        std::memmove(M0, M1, rowsize * sizeof(L0[0]));
        std::memmove(M1, M2, rowsize * sizeof(L0[0]));
        h += 1;
    }
    return true;
}

const std::string five_prime = "ATGGAGAGCCTTGTCCCTGG";
const std::string three_prime = "TTAACTTTAATCTCACATAG";
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

bool filter(const std::string& s) {
    if (s.length() < MIN_LENGTH) return false;
    if (s.length() > MAX_LENGTH) return false;
    size_t n = std::count(s.begin(), s.end(), 'N');
    float perc_gaps = n/ (float) s.length();
    if (perc_gaps > 0.2) {
        std::cerr << "too many gaps'" << std::endl; return false;
    }
    return true;
}

void DistCalculator::init_state_array(int* state_arr, int rowsize) {
    for (int i = 0; i < rowsize; ++i) {
        state_arr[i] = -2;
    }
    for (int i = rowsize; i < 2*rowsize; ++i) {
        state_arr[i] = -1;
    }
    for (int i = 2*rowsize; i < 3*rowsize; ++i) {
        state_arr[i] = -1;
    }
    for (int i = 3*rowsize; i < 4*rowsize; ++i) {
        state_arr[i] = 0;
    }
    for (int i = 4*rowsize; i < 5*rowsize; ++i) {
        state_arr[i] = 0;
    }
    for (int i = 5*rowsize; i < 6*rowsize; ++i) {
        state_arr[i] = 0;
    }
}

void DistCalculator::init_state_triple(int* state_triple, int len1, int len2) {
    state_triple[0] = 0;
    state_triple[1] = -2 * len1;
    state_triple[2] = 2 * len2;
}

void DistCalculator::query_samples_against_refs(std::string sample_fasta_fname, std::string ref_fasta_fname) {
    std::vector<std::pair<std::string, std::string>> queries = read_fasta(sample_fasta_fname);
    std::vector<std::pair<std::string, std::string>> refs = read_fasta(ref_fasta_fname);   
    int state_arr[6*MAX_ROW_SIZE] = {0};
    int state_triple[3] = {0,0,0};
    for (auto& p1: queries) {
        for (auto& p2: refs) {
            int rowsize = p1.second.length() + p2.second.length() + 1;
            init_state_array(state_arr, rowsize);
            init_state_triple(state_triple, p1.second.length(), p2.second.length());
            int mm_ind = 5*rowsize + p2.second.length();
            calculate_dist(p1.second, p2.second, state_triple, state_arr, rowsize, 5, 100, false);
            std::cout << p1.first << "," << p2.first << "," << state_triple[0] << "," << state_arr[mm_ind] << std::endl;
        }
    }
}

std::vector<std::pair<std::string, std::string>> DistCalculator::read_fasta(std::string fasta_fname) {
    std::cerr << "Reading " << fasta_fname << std::endl;
    std::vector<std::pair<std::string,std::string>> records;
    std::string curr_str = "";
    std::string curr_h = "";
    std::string l = "";
    int rec_counter = 0;
    int c = 0;
    std::ifstream fasta_file(fasta_fname);
    while (std::getline(fasta_file,l)) {
        l.erase(std::remove(l.begin(), l.end(), '\n'), l.end());
        l.erase(std::remove(l.begin(), l.end(), '\r'), l.end());
        if (l[0] == '>') {
            if (c > 0) {
                std::transform(curr_str.begin(), curr_str.end(), curr_str.begin(),
               [](unsigned char c){ return std::toupper(c); });
                convert_nonstandard_to_N(curr_str);
                bool cds = trim_to_cds(curr_str);
                if (cds && filter(curr_str)) {
                    std::pair<std::string,std::string> record(curr_h, curr_str);
                    records.push_back(record);
                }
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
    std::transform(curr_str.begin(), curr_str.end(), curr_str.begin(),
           [](unsigned char c){ return std::toupper(c); });
    convert_nonstandard_to_N(curr_str);
    bool cds = trim_to_cds(curr_str);
    if (cds && filter(curr_str)) {
        std::pair<std::string,std::string> record(curr_h, curr_str);
        records.push_back(record);
    }
    std::cerr << "Excluded " << rec_counter-records.size() << " of " << records.size() << " records due to filter..." << std::endl;
    return records;
}


