#include "uk2.hpp"
#include <string>


int main(int argc, char* argv[]) {
    int snpmax = atoi(argv[3]);
    DistCalculator dc = DistCalculator();
    dc.query_samples_against_refs(argv[1], argv[2], snpmax);
}

