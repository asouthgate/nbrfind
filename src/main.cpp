#include "uk2.hpp"
#include <string>


int main(int argc, char* argv[]) {
    int snpmax = atoi(argv[3]);
    int max_slide = atoi(argv[4]);
    DistCalculator dc = DistCalculator();
    dc.query_samples_against_refs(argv[1], argv[2], snpmax, max_slide);
}

