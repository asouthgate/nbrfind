StateData(int* input_state_quintuple, int* input_state_arr, int rowsize) {
    state_arr = input_state_arr;
    state_quintuple = input_state_quintuple;
    
    lower_bound = state_quintuple[1];
    upper_bound = state_quintuple[2];
    // the diagonal to resume on, at (h,d)
    d_resume = state_quintuple[3];
    // What is maxi_resume?
    maxi_resume = state_quintuple[4];
    if (d_resume > 0) assert(maxi_resume > 0);
    // Init L arrays
    L0 = state_arr;
    L1 = state_arr + rowsize;
    L2 = state_arr + 2 * rowsize;
    // Init M arrays; number of mismatch operations
    M0 = state_arr + 3 * rowsize;
    M1 = state_arr + 4 * rowsize;
    M2 = state_arr + 5 * rowsize;
    // Auxiliary arrays:
    // Init NM arrays; number of match operations
    NM0 = state_arr + 6 * rowsize;
    NM1 = state_arr + 7 * rowsize;
    NM2 = state_arr + 8 * rowsize;
    // Init NN arrays; number of `N' match operations
    NN0 = state_arr + 9 * rowsize;
    NN1 = state_arr + 10 * rowsize;
    NN2 = state_arr + 11 * rowsize;
    // starting h
    h = state_quintuple[0];
}

void StateData::freeze(int prev_lower_bound, int prev_upper_bound, int d, int maxi) {
    lower_bound = prev_lower_bound;
    upper_bound = prev_upper_bound;
    d_resume = d;
    maxi_resume = maxi;
    std::memmove(L2, L1, rowsize * sizeof(L0[0]));
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
