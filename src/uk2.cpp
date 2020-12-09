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


MoveResult cal_move(StateData& sd, const int& d, const int& ld, const int& m, const int& n, const int& imax) {
    int lmove = sd.L0[ld-1];
    int matchmove = std::min(sd.L1[ld] + 1, imax);
    int rmove = std::min(sd.L0[ld+1] + 1, imax);
    int prev_mm = sd.M1[ld];
    int maxi = -1;
    int prev_NM, prev_NN;
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
    std::cerr << "maxi: " << maxi << std::endl;
    return MoveResult(maxi, prev_NM, prev_NN);
}

bool DistCalculator::calculate_dist_sd(std::string s1, std::string s2, StateData& sd, int snpmax, int slide_threshold, bool freeze) {
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
    while ((sd.h) < maxh) {
        sd.print_debug();
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
             dstart = std::max(-(sd.h/2), sd.lower_bound);
        }
        // Why is this dmax?
        int dmax = std::min( (sd.h/2), sd.upper_bound) + 1;
        for (int d = dstart; d < dmax; ++d) {
            std::cerr << d << std::endl;
            // Get diagonal index
            int ld = m + d;
            // Why? what is imax?
            // imax is the maximum i that a diagonal can get to
            int imax = imax_arr[ld];
            assert (ld >= 0);
            assert (sd.lower_bound + m >= 0);
            assert (sd.upper_bound + m <= m + n + 1);
            assert (sd.L1[ld] < imax);
            // Construct move result to hold maxi, prev_NM, prev_NN
            MoveResult move_result;

            // MOVE 1
            if (sd.maxi_resume == 0) {
                move_result = cal_move(sd, d, ld, m, n, imax);
            }
            else { 
                move_result = MoveResult(sd.maxi_resume, sd.NM2[ld], sd.NN2[ld]);
                sd.maxi_resume = 0;
            }
            
            // MOVE 2: slide
            int sl, mn; 
            if (move_result.maxi < imax) {
                std::pair<int,int> res =  slide(d, move_result.maxi, s1, s2, imax);
                sl = res.first; mn = res.second;
                assert (sl <= imax);
                move_result.prev_NM += (sl-move_result.maxi);
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
                    sd.freeze(prev_lower_bound, prev_upper_bound, d, move_result.maxi);
                    return true;
                }
            }

            // UPDATE BOUNDS
            assert (sd.M2[ld] >= 0);
            // Update bounds if hit imax
            if (ld <= n && sd.L2[ld] >= imax) {
                sd.lower_bound = d + 1;
                assert (sd.lower_bound + m >= 0);
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
        sd.print_debug();

        // SWAP POINTERS
        sd.swap_pointers();
        sd.h += 1;
    }
    return true;
}

void DistCalculator::query_samples_against_refs(std::string sample_fasta_fname, std::string ref_fasta_fname, int k, double epsilon) {
    std::vector<std::pair<std::string, std::string>> queries = read_fasta(sample_fasta_fname, MIN_LENGTH, MAX_LENGTH);
    std::vector<std::pair<std::string, std::string>> refs = read_fasta(ref_fasta_fname, MIN_LENGTH, MAX_LENGTH);
    size_t max_seq_len = 0;
    for (int si = 0; si < refs.size(); ++si) {
        max_seq_len = std::max(max_seq_len, refs[si].second.length());
    }
    for (int si = 0; si < refs.size(); ++si) {
        max_seq_len = std::max(max_seq_len, queries[si].second.length());
    }
    StateData sd(max_seq_len*2 + 1);
    for (auto& p1 : queries) {
        for (auto& p2 : refs) {
            int rowsize = p1.second.length() + p2.second.length() + 1;
            sd.init_state_quintuple(p1.second.length(), p2.second.length());
            sd.init_state_array(rowsize);
            int n = p2.second.length();
            bool res = calculate_dist_sd(p1.second, p2.second, sd, k, 100, false);
//            double p = betap(state_arr[mm_ind], state_arr[NM_ind], state_arr[NN_ind], k);
            std::cout << p1.first << "," << p2.first << "," << sd.h << "," << sd.M2[n] << ","
                          << sd.NM2[n] << "," << sd.NN2[n] << std::endl;
        }
    }
}


