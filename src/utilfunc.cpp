#include <algorithm> 
#include <string.h>
#include <iostream>

extern "C" {
    int max_match_ind(const char* a, const char* b) {
//        std::cout << strlen(a) << " " << strlen(b) << std::endl;
        for (int j = 0; j < std::min(strlen(a), strlen(b)); ++j) {
//            std::cout << a[j] << " " << b[j] << std::endl;
            if (a[j] != b[j]) return j-1;
        }
        return strlen(a)-1;
    }
}
