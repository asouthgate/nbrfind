#include "sdata.hpp"
#include <assert.h>
#include <cstring>

StateData::StateData(int MAX_ROW_SIZE_) {
    MAX_ROW_SIZE = MAX_ROW_SIZE_;
    state_arr = new int[12*MAX_ROW_SIZE]();
    
    if (d_resume > 0) assert(maxi_resume > 0);
    // Init L arrays
    L0 = state_arr;
    L1 = state_arr + MAX_ROW_SIZE;
    L2 = state_arr + 2 * MAX_ROW_SIZE;
    // Init M arrays; number of mismatch operations
    M0 = state_arr + 3 * MAX_ROW_SIZE;
    M1 = state_arr + 4 * MAX_ROW_SIZE;
    M2 = state_arr + 5 * MAX_ROW_SIZE;
    // Auxiliary arrays:
    // Init NM arrays; number of match operations
    NM0 = state_arr + 6 * MAX_ROW_SIZE;
    NM1 = state_arr + 7 * MAX_ROW_SIZE;
    NM2 = state_arr + 8 * MAX_ROW_SIZE;
    // Init NN arrays; number of `N' match operations
    NN0 = state_arr + 9 * MAX_ROW_SIZE;
    NN1 = state_arr + 10 * MAX_ROW_SIZE;
    NN2 = state_arr + 11 * MAX_ROW_SIZE;
}

void StateData::freeze(int prev_lower_bound, int prev_upper_bound, int d, int maxi) {
    lower_bound = prev_lower_bound;
    upper_bound = prev_upper_bound;
    d_resume = d;
    maxi_resume = maxi;
    std::memmove(L2, L1, MAX_ROW_SIZE * sizeof(L0[0]));
}

void StateData::swap_pointers() {
    L0 = L1;
    L1 = L2;
    M0 = M1;
    M1 = M2;
    NM0 = NM1;
    NM1 = NM2;
    NN0 = NN1;
    NN1 = NN2;
}

void StateData::init_state_array(int rowsize) {
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

void StateData::init_state_quintuple(int len1, int len2) {
    // starting h
    h = 0;
    lower_bound = -len1;
    upper_bound = len2;
    // the diagonal to resume on, at (h,d)
    d_resume = 0;
    // What is maxi_resume?
    maxi_resume = 0;
}
