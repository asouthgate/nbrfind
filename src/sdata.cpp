#include "sdata.hpp"
#include <assert.h>
#include <cstring>
#include <iostream>

void StateData::print_debug() {
    std::cerr << "h: " << h << std::endl;
    print_arrays();
    std::cerr << "lb: " << lower_bound << "ub: " << upper_bound << std::endl;
    std::cerr << std::endl;

}

void StateData::print_arrays() {
    for (int i = 0; i < MAX_ROW_SIZE; ++i) {
        std::cerr << L0[i] << " ";
    }
    std::cerr << std::endl;
    for (int i = 0; i < MAX_ROW_SIZE; ++i) {
        std::cerr << L1[i] << " ";
    }
    std::cerr << std::endl;
    for (int i = 0; i < MAX_ROW_SIZE; ++i) {
        std::cerr << L2[i] << " ";
    }
    std::cerr << std::endl;
    for (int i = 0; i < MAX_ROW_SIZE; ++i) {
        std::cerr << M0[i] << " ";
    }
    std::cerr << std::endl;
    for (int i = 0; i < MAX_ROW_SIZE; ++i) {
        std::cerr << M1[i] << " ";
    }
    std::cerr << std::endl;
    for (int i = 0; i < MAX_ROW_SIZE; ++i) {
        std::cerr << M2[i] << " ";
    }
    std::cerr << std::endl;
}

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
    // 0 <- 1
    // 1 <- 2
    // 2 <- 0 AND ENSURE NO MINUS 2S PRESENT IF h=0
    int* tmp = L0;
    L0 = L1;
    L1 = L2;
    L2 = tmp;
    tmp = M0;
    M0 = M1;
    M1 = M2;
    M2 = tmp;
    tmp = NM0;
    NM0 = NM1;
    NM1 = NM2;
    NM2 = tmp;
    tmp = NN0;
    NN0 = NN1;
    NN1 = NN2;
    NN2 = tmp;
//    if (h == 0) {
//        for (int i = 0; i < MAX_ROW_SIZE; ++i) {
//            L2[i] = -1;
//        }
//    }
}

//void StateData::init_state_array(int rowsize) {
//    for (int i = 0; i < rowsize; ++i) {
//        state_arr[i] = -2;
//    }
//    for (int i = rowsize; i < 3*rowsize; ++i) {
//        state_arr[i] = -1;
//    }
//    for (int i = 3*rowsize; i < 12*rowsize; ++i) {
//        state_arr[i] = 0;
//    }
//
//}

void StateData::init_state_array(int rowsize) {
    // Rowsize is the current amount of each array we will need to fill
    for (int i = 0; i < rowsize; ++i) {
        L0[i] = -1;
        L1[i] = -1; L2[i] = -1;
        M0[i] = 0; M1[i] = 0; M2[i] = 0;
        NM0[i] = 0; NM1[i] = 0; NM2[i] = 0;
        NN0[i] = 0; NN1[i] = 0; NN2[i] = 0;
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
