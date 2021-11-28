#include <algorithm> 
#include <string.h>
#include <iostream>

extern "C" {
/**
* Find the maximum matching index between two C strings.
*
* @param first the first string.
* @param second the second string.
* @return the maximum index that the strings match with.
*/
    int max_match_ind(const char* first, const char* second) {
        for (int j = 0; j < std::min(strlen(first), strlen(second)); ++j) {
            if (first[j] != second[j]) return j-1;
        }
        return strlen(a)-1;
    }
}

void cerrarr(int* v, int rowsize) {
    for (int pj = 0; pj < rowsize; ++pj) {
        std::cerr << v[pj] << " ";
    }
    std::cerr << std::endl;
}

