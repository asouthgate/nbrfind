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

std::pair<int,int> slide(int d, int i, const std::string& s1, const std::string& s2, int imax) {
    assert (i < imax);
    int Nc = 0;
    int qend = std::min(s1.length(), s2.length()-d);
    // Fastest way to optimize this?
    for (int q = i; q < qend; ++q) {
        if ((s1[q] != 'N') && (s2[q+d] != 'N')){ 
            if ((s1[q] != s2[q+d])) return std::pair<int,int>(q,Nc);
        }
        else {
            Nc += 1;
        }
    }
    return std::pair<int,int>(imax, Nc);
}

void cerrarr(int* v, int rowsize) {
    for (int pj = 0; pj < rowsize; ++pj) {
        std::cerr << v[pj] << " ";
    }
    std::cerr << std::endl;
}

unsigned long long binom(unsigned int n, unsigned int k) {
    if (n == k || k == 0) return 1;
    return binom(n-1, k-1) * n/k;
}

unsigned long long beta(unsigned long long u, unsigned long long v) { 
    return std::tgammal(u)*std::tgammal(v)/std::tgamma(u+v);
}

double betap(int d, int M, int N, int k, int a=1, int b=1) {
    if (d >= k) return 0;
    double sum = 0;
    for (int z = 0; z < k; ++z) {
        unsigned long long b = binom(N,z);
        long double numo = beta(z-a+d,N-z+b+M-d);
        long double deno = beta(a+d, b+M-d);
    }
} 

int cal_imax(int d, int m, int n) {
    if (m >= n) return d <= 0 ? m : m-d;
    int v = n-m;
    return d <= v ? m : m-d+v;
}

std::vector<int> cal_imax_arr(int m, int n) {
    std::vector<int> res;
    for (int ld = 0; ld < m+n+1; ++ld) {
        res.push_back(cal_imax(ld-m,m,n));
    }
    return res;
}

/*
* Compute the diagonal dynamic programming matrix and auxiliaries
*
*
*/
bool DistCalculator::calculate_dist(std::string s1, std::string s2, int* state_triple, int* state_arr, int rowsize, int snpmax, int slide_threshold, bool freeze) {
    int m = s1.length();
    int n = s2.length();
    int h_start = state_triple[0];
    int lower_bound = state_triple[1];
    int upper_bound = state_triple[2];
    // Init L arrays
    auto L0 = state_arr;
    auto L1 = state_arr + rowsize;
    auto L2 = state_arr + 2*rowsize;
    // Init M arrays; number of mismatch operations
    auto M0 = state_arr + 3*rowsize;
    auto M1 = state_arr + 4*rowsize;
    auto M2 = state_arr + 5*rowsize;
    // Auxiliary arrays:
    // Init NM arrays; number of match operations
    auto NM0 = state_arr + 6*rowsize;
    auto NM1 = state_arr + 7*rowsize;
    auto NM2 = state_arr + 8*rowsize;
    // Init NN arrays; number of `N' match operations
    auto NN0 = state_arr + 9*rowsize;
    auto NN1 = state_arr + 10*rowsize;
    auto NN2 = state_arr + 11*rowsize;

    if (freeze && n == 0) { return false; }

    int h = h_start;
    std::vector<int> imax_arr = cal_imax_arr(m, n);
    while (h < 2*(m+n)+1) {
//        cerrarr(NM0,rowsize); cerrarr(NM1,rowsize); cerrarr(NM2,rowsize); 
//        std::cerr << h << std::endl;
        state_triple[0] = h;
        int prev_lower_bound = lower_bound;
        int prev_upper_bound = upper_bound;
        for (int d = std::max(-(h/2), lower_bound); d < std::min((h/2),upper_bound)+1; ++d) {
            int ld = m+d;
            int imax = imax_arr[ld];
            assert (ld >= 0);
            assert (lower_bound + m >= 0);
            assert (upper_bound + m <= m+n+1);
            assert (L1[ld] < imax);
            int maxi = -1;
            int lmove = L0[ld-1];
            int matchmove = std::min(L1[ld] + 1, imax);
            int rmove = std::min(L0[ld+1] + 1, imax);
            int prev_mm = M1[ld];
            int prev_NM = -1;
            int prev_NN = -1;
            // TODO: check the bounds are not redundant
            if ((d > -h) && (ld-1 > 0)) {
                maxi = std::max(maxi, lmove);
                M2[ld] = M0[ld-1];
                prev_NM = NM0[ld-1];
                prev_NN = NN0[ld-1];
            }
            if ((d < h) && (ld+1 < m+n+1)) {
                maxi = std::max(maxi, rmove);
                if (rmove >= lmove) { 
                    M2[ld] = M0[ld+1]; 
                    prev_NM = NM0[ld+1];
                    prev_NN = NN0[ld+1];
                }
            }
            if (matchmove >= maxi) {
                maxi = matchmove;
                prev_NM = NM1[ld];
                prev_NN = NN1[ld];
                if (maxi > 0) {
                    prev_NM += 1;
                    M2[ld] = prev_mm + 1;
                }
            }
            // TODO: redundant?
            maxi = std::min(maxi, imax);
            int sl, mn;
            if (maxi < imax) {
                std::pair<int,int> res =  slide(d, maxi, s1, s2, imax);
                sl = res.first; mn = res.second;
                assert (sl <= imax);
                prev_NM += (sl-maxi);
                prev_NN += mn;
                L2[ld] = sl;
            }
            else { L2[ld] = maxi; }
            NM2[ld] = prev_NM;
            NN2[ld] = prev_NN;
            assert (maxi <= imax);
            assert (L2[ld] <= imax);
                    
            if (freeze && d >= n-m) {
                if (L2[ld] == imax) {
                    return true;
                }
            }

            assert (M2[ld] >= 0);
            if (ld <= n && L2[ld] >= imax) {
                lower_bound = d+1; state_triple[1] = d;
                assert (lower_bound + m >= 0);
            }
            if (ld > n && L2[ld] >= imax) {
                upper_bound = d-1; state_triple[2] = d-1;
                assert (upper_bound + m >= n);
                assert (upper_bound + m <= n+m+1);
                break;
            } 
            if (L2[n] == m) {
                return true;
            }
        }
        // TODO: memcopy? don't overlap in range
        std::memmove(L0, L1, rowsize * sizeof(L0[0]));
        std::memmove(L1, L2, rowsize * sizeof(L0[0]));

        std::memmove(M0, M1, rowsize * sizeof(L0[0]));
        std::memmove(M1, M2, rowsize * sizeof(L0[0]));

        std::memmove(NM0, NM1, rowsize * sizeof(L0[0]));
        std::memmove(NM1, NM2, rowsize * sizeof(L0[0]));

        std::memmove(NN0, NN1, rowsize * sizeof(L0[0]));
        std::memmove(NN1, NN2, rowsize * sizeof(L0[0]));
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
    int state_arr[12*MAX_ROW_SIZE] = {0};
    int state_triple[3] = {0,0,0};
    for (auto& p1: queries) {
        for (auto& p2: refs) {
            int rowsize = p1.second.length() + p2.second.length() + 1;
            init_state_array(state_arr, rowsize);
            init_state_triple(state_triple, p1.second.length(), p2.second.length());
            int mm_ind = 5*rowsize + p2.second.length();
            int NM_ind = 8*rowsize + p2.second.length();
            int NN_ind = 11*rowsize + p2.second.length();
            calculate_dist(p1.second, p2.second, state_triple, state_arr, rowsize, 5, 100, false);
            std::cout << p1.first << "," << p2.first << "," 
                      << state_triple[0] << "," << state_arr[mm_ind] 
                      << state_arr[NM_ind] << "," << state_arr[NN_ind] << std::endl;
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


