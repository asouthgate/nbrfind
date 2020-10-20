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
#include "sdata.hpp"

std::pair<int, int> slide(int d, int i, const std::string& s1, const std::string& s2, int imax) {
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


MoveResult cal_move(StateData& state_data, const int& d, const int& ld, const int& m, const int& n, const int& imax) {
    int lmove = sd.L0[ld-1];
    int matchmove = std::min(sd.L1[ld] + 1, imax);
    int rmove = std::min(sd.L0[ld+1] + 1, imax);
    int prev_mm = sd.M1[ld];
    // TODO: check the bounds are not redundant
    if ((d > -sd.h) && (ld - 1 > 0)) {
        maxi = std::max(maxi, lmove);
        sd.M2[ld] = sd.M0[ld-1];
        prev_NM = sd.NM0[ld-1];
        prev_NN = sd.NN0[ld-1];
    }
    if ((d < sd.h) && (ld + 1 < m + n + 1)) {
        maxi = std::max(maxi, rmove);
        if (rmove >= lmove) { 
            sd.M2[ld] = sd.M0[ld+1]; 
            prev_NM = sd.NM0[ld+1];
            prev_NN = sd.NN0[ld+1];
        }
    }
    if (matchmove >= maxi) {
        maxi = matchmove;
        prev_NM = sd.NM1[ld];
        prev_NN = sd.NN1[ld];
        if (maxi > 0) {
            prev_NM += 1;
            sd.M2[ld] = prev_mm + 1;
        }
    }
    // TODO: redundant?
    maxi = std::min(maxi, imax);
    return MoveResult(maxi, prev_NM, prev_NN);
}

bool DistCalculator::calculate_dist_sd(std::string s1, std::string s2, StateData& state_data, int snpmax, int slide_threshold, bool freeze) {
    // Initialization
    int m = s1.length();
    int n = s2.length();
    if (freeze && n == 0) { return true; }
    // What is the imax_arr for?
    std::vector<int> imax_arr;
    imax_arr.reserve(m + n + 1);
    fill_imax_arr(imax_arr, m, n);
    // Why is this maxh?
    int maxh = (2 * (m + n) + 1);
    while ( (sd.h) < maxh) {
        // Why do we store prev_lower and prev_upper?
        int prev_lower_bound = sd.lower_bound;
        int prev_upper_bound = sd.upper_bound;
        int dstart;
        // Find first diagonal to start on
        if (sd.d_resume > 0) { 
            dstart = sd.d_resume;
            sd.d_resume = 0;
        }
        else {
             dstart = std::max(-(h/2), lower_bound);
        }
        // Why is this dmax?
        int dmax = std::min( (sd.h/2), sd.upper_bound) + 1;
        for (int d = dstart; d < dmax; ++d) {
            // Get diagonal index
            int ld = m + d;
            // Why? what is imax?
            int imax = imax_arr[ld];
            assert (ld >= 0);
            assert (sd.lower_bound + m >= 0);
            assert (sd.upper_bound + m <= m + n + 1);
            assert (sd.L1[ld] < imax);
            // Construct move result to hold maxi, prev_NM, prev_NN
            MoveResult move_result = move_result();

            // MOVE 1
            if (sd.maxi_resume == 0) {
                cal_move(move_result);
            }
            else { 
                move_result = MoveResult(sd.maxi_resume, sd.NM2[ld], sd.NN2[ld]);
                sd.maxi_resume = 0;
            }
            
            // MOVE 2: slide
            int sl, mn; 
            // What is maxi vs imax??
            if (move_result.maxi < imax) {
                std::pair<int,int> res =  slide(d, move_result.maxi, s1, s2, imax);
                sl = res.first; mn = res.second;
                assert (sl <= imax);
                move_result.prev_NM += (sl-maxi);
                // heuristic abandonment
                if (sl - move_result.maxi > slide_threshold) {
                    if (sd.M2[ld] >= snpmax) { return false; }
                }
                move_result.prev_NN += mn;
                sd.L2[ld] = sl;
            }
            else { sd.L2[ld] = move_result.maxi; }

            sd.NM2[ld] = move_result.prev_NM;
            sd.NN2[ld] = move_result.prev_NN;
            assert(move_result.maxi <= imax);
            assert(sd.L2[ld] <= imax);
            //****
                    
            // FREEZE CHECK
            if (freeze && d >= n-m) {
                if (sd.L2[ld] == imax) {
                    sd.freeze(prev_lower_bound, prev_upper_bound, d, maxi);
                    return true;
                }
            }

            // UPDATE BOUNDS
            assert (sd.M2[ld] >= 0);
            // Update bounds if hit imax
            if (ld <= n && sd.L2[ld] >= imax) {
                sd.lower_bound = d + 1;
                assert (lower_bound + m >= 0);
            }
            if (ld > n && sd.L2[ld] >= imax) {
                sd.upper_bound = d - 1;
                assert (sd.upper_bound + m >= n);
                assert (sd.upper_bound + m <= n + m + 1);
                break;
            } 
            if (sd.L2[n] == m) {
                return true;
            }
        }

        // SWAP POINTERS
        sd.swap_pointers();
        sd.h += 1;
    }
    return true;
}



/*
* Compute the diagonal dynamic programming matrix and auxiliaries.
*
* This program works in place on the inputted state.
*/
bool DistCalculator::calculate_dist(std::string s1, std::string s2, int* state_quintuple, int* state_arr, int rowsize, int snpmax, int slide_threshold, bool freeze) {

    // load_state()
    // Initialization and creating references into data; 
    int m = s1.length();
    int n = s2.length();
    if (freeze && n == 0) { return true; }
    int lower_bound = state_quintuple[1];
    int upper_bound = state_quintuple[2];
    // the diagonal to resume on, at (h,d)
    int d_resume = state_quintuple[3];
    int maxi_resume = state_quintuple[4];
    if (d_resume > 0) assert(maxi_resume > 0);
    // Init L arrays
    auto L0 = state_arr;
    auto L1 = state_arr + rowsize;
    auto L2 = state_arr + 2 * rowsize;
    // Init M arrays; number of mismatch operations
    auto M0 = state_arr + 3 * rowsize;
    auto M1 = state_arr + 4 * rowsize;
    auto M2 = state_arr + 5 * rowsize;
    // Auxiliary arrays:
    // Init NM arrays; number of match operations
    auto NM0 = state_arr + 6 * rowsize;
    auto NM1 = state_arr + 7 * rowsize;
    auto NM2 = state_arr + 8 * rowsize;
    // Init NN arrays; number of `N' match operations
    auto NN0 = state_arr + 9 * rowsize;
    auto NN1 = state_arr + 10 * rowsize;
    auto NN2 = state_arr + 11 * rowsize;
    // starting h
    int h = state_quintuple[0];

    // What is the imax_arr for?
    std::vector<int> imax_arr;
    imax_arr.reserve(m + n + 1);
    fill_imax_arr(imax_arr, m, n);
    // Move calc to maxh
    while ( h < (2 * (m + n) + 1) ) {
        // Update the state
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

        // Swap the pointers! We only hold a few in memory
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


