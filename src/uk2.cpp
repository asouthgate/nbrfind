#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1
#include <cmath>
#include <assert.h>
#include <math.h>
#include <algorithm>
#include <string>
#include <iostream>
#include <cstring>
#include <fstream>
#include <set>
#include "uk2.hpp"
#include "probfunc.hpp"
#include "seqio.hpp"

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

int cal_imax(int d, int m, int n) {
    if (m >= n) return d <= 0 ? m : m-d;
    int v = n-m;
    return d <= v ? m : m-d+v;
}

void fill_imax_arr(std::vector<int> & v, int m, int n) {
    for (int ld = 0; ld < m+n+1; ++ld) {
        v[ld] = cal_imax(ld-m,m,n);
    }
}

/*
* Compute the diagonal dynamic programming matrix and auxiliaries
*
*
*/
bool DistCalculator::calculate_dist(std::string s1, std::string s2, int* state_quintuple, int* state_arr, int rowsize, int snpmax, int slide_threshold, bool freeze) {
    int m = s1.length();
    int n = s2.length();
    int h_start = state_quintuple[0];
    int lower_bound = state_quintuple[1];
    int upper_bound = state_quintuple[2];
    // the diagonal to resume on, at (h,d)
    int d_resume = state_quintuple[3];
    int maxi_resume = state_quintuple[4];
    if (d_resume > 0) assert(maxi_resume > 0);
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

    if (freeze && n == 0) { return true; }

    int h = h_start;
    std::vector<int> imax_arr;
    imax_arr.reserve(m+n+1);
    fill_imax_arr(imax_arr, m, n);
    while (h < 2*(m+n)+1) {
//        cerrarr(L0,rowsize); cerrarr(L1,rowsize); cerrarr(L2,rowsize); 
//        std::cerr << h << std::endl;
        state_quintuple[0] = h;
        int prev_lower_bound = lower_bound;
        int prev_upper_bound = upper_bound;
        int dstart;
        if (d_resume > 0) { 
            dstart = d_resume;
            d_resume = 0;
        }
        else {
             dstart = std::max(-(h/2), lower_bound);
        }
        for (int d = dstart; d < std::min((h/2),upper_bound)+1; ++d) {
            int ld = m+d;
            int imax = imax_arr[ld];
            assert (ld >= 0);
            assert (lower_bound + m >= 0);
            assert (upper_bound + m <= m+n+1);
            assert (L1[ld] < imax);
            int maxi = -1;
            int prev_NM = -1;
            int prev_NN = -1;
            if (maxi_resume == 0) {
                int lmove = L0[ld-1];
                int matchmove = std::min(L1[ld] + 1, imax);
                int rmove = std::min(L0[ld+1] + 1, imax);
                int prev_mm = M1[ld];
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
            }
            else { 
                maxi = maxi_resume;
                maxi_resume = 0;
                prev_NM = NM2[ld];
                prev_NN = NN2[ld];
            }
            int sl, mn;
            if (maxi < imax) {
                std::pair<int,int> res =  slide(d, maxi, s1, s2, imax);
                sl = res.first; mn = res.second;
                assert (sl <= imax);
                prev_NM += (sl-maxi);
                // heuristic abandonment
                if (sl-maxi > slide_threshold) {
                    if (M2[ld] >= snpmax) { return false; }
                }
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
//                    std::cerr << "freezing on" << std::endl;
//                    cerrarr(L0,rowsize); cerrarr(L1,rowsize); cerrarr(L2,rowsize); 
                    state_quintuple[1] = prev_lower_bound;
                    state_quintuple[2] = prev_upper_bound;
                    state_quintuple[3] = d;
                    state_quintuple[4] = maxi;
                    std::memmove(L2, L1, rowsize * sizeof(L0[0]));
                    return true;
                }
            }

            assert (M2[ld] >= 0);
            if (ld <= n && L2[ld] >= imax) {
                lower_bound = d+1; state_quintuple[1] = d+1;
                assert (lower_bound + m >= 0);
            }
            if (ld > n && L2[ld] >= imax) {
                upper_bound = d-1; state_quintuple[2] = d-1;
                assert (upper_bound + m >= n);
                assert (upper_bound + m <= n+m+1);
                break;
            } 
            if (L2[n] == m) {
                return true;
            }
        }

        // Don't need to memcopy! Swap pointers;
        L0 = L1;
        L1 = L2;
        M0 = M1;
        M1 = M2;
        NM0 = NM1;
        NM1 = NM2;
        NN0 = NN1;
        NN1 = NN2;
        h += 1;
    }
    return true;
}

void DistCalculator::init_state_array(int* state_arr, int rowsize) {
    for (int i = 0; i < rowsize; ++i) {
        state_arr[i] = -2;
    }
    for (int i = rowsize; i < 3*rowsize; ++i) {
        state_arr[i] = -1;
    }
    for (int i = 3*rowsize; i < 12*rowsize; ++i) {
        state_arr[i] = 0;
    }

}

void DistCalculator::init_state_quintuple(int* state_quintuple, int len1, int len2) {
    state_quintuple[0] = 0;
    state_quintuple[1] = -len1;
    state_quintuple[2] = len2;
    state_quintuple[3] = 0;
    state_quintuple[4] = 0;
}

bool pair_compare(const std::pair<std::string,std::string> &a, const std::pair<std::string,std::string> &b) {
    return a.second < b.second;
}

int get_last_prefix_match(const std::string & s1, const std::string & s2) {
    for (int j = 0; j < std::min(s1.length(), s2.length()); ++j) {
         if (s1[j] != s2[j]) return j-1;
    }
    return std::min(s1.length(), s2.length())-1;
}

std::vector<int> get_prefix_array(const std::vector<std::pair<std::string,std::string>>& records) {
    std::vector<int> res;
    for (int j = 1; j < records.size(); ++j) {
        res.push_back(get_last_prefix_match(records[j].second, records[j-1].second));
    }
    return res;
}

std::vector<Cluster> DistCalculator::get_clusters(const std::vector<std::pair<std::string,std::string>>& tmprefs, int k, double epsilon) {
    std::cerr << "Fetching clusters" << std::endl;
    std::vector<Cluster> clusters;
    // Initiate available set
    std::vector<int> available_members;
    // init required arrays for comparisons
    int* state_arr = new int[12*MAX_ROW_SIZE]();
    int state_quintuple[5] = {0,0,0,0,0};
    for (int i = 1; i < tmprefs.size(); ++i) available_members.push_back(i);
    while (available_members.size() > 0) {
        Cluster cluster;
        // take a new curr_centroid
        cluster.xi = available_members[0];
        for (int j = 1; j < available_members.size(); ++j) {
            int refi = available_members[j];
            std::pair<std::string, std::string> p1 = tmprefs[refi];
            std::pair<std::string, std::string> p2 = tmprefs[cluster.xi];
            int rowsize = p1.second.length() + p2.second.length() + 1;
            init_state_quintuple(state_quintuple, p1.second.length(), p2.second.length());
            init_state_array(state_arr, rowsize);
            int mm_ind = 5*rowsize + p2.second.length();
            int NM_ind = 8*rowsize + p2.second.length();
            int NN_ind = 11*rowsize + p2.second.length();
            bool res = calculate_dist(p1.second, p2.second, state_quintuple, state_arr, rowsize, k, 100, false);
            double p = betap(state_arr[mm_ind], state_arr[NM_ind], state_arr[NN_ind], 2*k);
            if (res && p > epsilon) {
                cluster.members.push_back(refi);
            }
        }
        std::cerr << "creating new cluster of size "<< cluster.members.size() << std::endl;
        // Now delete the cluster members from available members; go backwards
        std::cerr << "deleting... " << " ";
        available_members.erase(std::remove(available_members.begin(), available_members.end(), cluster.xi), available_members.end());
        for (int cmj : cluster.members) {
            std::cerr << cmj << " ";
            available_members.erase(std::remove(available_members.begin(), available_members.end(), cmj), available_members.end());
        }        
        std::cerr << std::endl;
        std::cerr << "remaining members "<< available_members.size() << std::endl;
        clusters.push_back(cluster);
    }
    delete[] state_arr;
    std::cerr << "Got " << clusters.size() << " clusters" << std::endl;
    return clusters;
}

void DistCalculator::query_samples_against_refs(std::string sample_fasta_fname, std::string ref_fasta_fname, int k, double epsilon) {
    std::vector<std::pair<std::string, std::string>> queries = read_fasta(sample_fasta_fname, MIN_LENGTH, MAX_LENGTH);
    std::vector<std::pair<std::string, std::string>> refs = read_fasta(ref_fasta_fname, MIN_LENGTH, MAX_LENGTH);
//    // Firstly create clusters and reps
//    std::cerr << "fetching clusters within " << k << " " << epsilon << std::endl;
//    std::vector<Cluster> clusters = get_clusters(refs, k, epsilon);
    // allocate these on the heap
    int* state_arr = new int[12*MAX_ROW_SIZE]();
    int state_quintuple[5] = {0,0,0,0,0};
    for (auto& p1 : queries) {
        for (auto& p2 : refs) {
//        for (Cluster& c : clusters) {
//            int xi = c.xi;
//            std::pair<std::string, std::string> p2 = refs[xi];
            int rowsize = p1.second.length() + p2.second.length() + 1;
            init_state_quintuple(state_quintuple, p1.second.length(), p2.second.length());
            init_state_array(state_arr, rowsize);
            int mm_ind = 5*rowsize + p2.second.length();
            int NM_ind = 8*rowsize + p2.second.length();
            int NN_ind = 11*rowsize + p2.second.length();
            bool res = calculate_dist(p1.second, p2.second, state_quintuple, state_arr, rowsize, k, 100, false);
            double p = betap(state_arr[mm_ind], state_arr[NM_ind], state_arr[NN_ind], k);
            std::cout << p1.first << "," << p2.first << "," << state_quintuple[0] << "," << state_arr[mm_ind] << ","
                          << state_arr[NM_ind] << "," << state_arr[NN_ind] << "," << p << std::endl;
        }
    }
    delete[] state_arr;
}


