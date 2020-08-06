#include "uk2.hpp"
#include <iostream>
#include <string.h>
#include <vector>
#include <cstddef>

extern "C" {
    void calculate_dist_resume_freeze(const char* a, const char* b, int* state_triple, int* state_arr, int rowsize) {
        DistCalculator dc = DistCalculator();
        dc.calculate_dist(a, b, state_triple, state_arr, rowsize, true);
    }

    void calculate_dist_resume_nofreeze(const char* a, const char* b, int* state_triple, int* state_arr, int rowsize) {
        DistCalculator dc = DistCalculator();
        dc.calculate_dist(a, b, state_triple, state_arr, rowsize, false);
    }
}

