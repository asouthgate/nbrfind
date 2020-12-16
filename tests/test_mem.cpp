#include<iostream>
#include<vector>
#include<string>
#include<assert.h>
#include"uk2.hpp"
#include"seqio.hpp"

int main(int argc, char* argv[]) {
    std::cerr << "testing" << std::endl; 
    std::string S = "NATATTNNATANANAAAATTTTTN";
    std::pair<std::vector<int>, std::string> res = strip_Ns(S);
    for (int j = 0; j < res.first.size(); ++j) {
        assert(res.second[j] == S[res.first[j]]);
    }
}
