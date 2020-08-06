#include "uk2.hpp"


int main(int argc, char* argv[]) {
    DistCalculator dc = DistCalculator();
    dc.query_samples_against_refs(argv[1], argv[2]);
}

