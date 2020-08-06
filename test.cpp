#include"uk2.hpp"
#include<iostream>

int main(int argc, char* argv[]) {
//    std::cout << argv[1] << " " << argv[2] << std::endl;
    Dist_calculator dc = Dist_calculator();
    std::vector<int> state = dc.calculate_dist(argv[1], argv[2], true);
    std::vector<int> final_state = dc.calculate_dist(argv[1], std::string(argv[2])+std::string(argv[3]), false, state);
    std::vector<int> final_state2 = dc.calculate_dist(argv[1], std::string(argv[2])+std::string(argv[3]), false, state);
    std::cout << final_state[0] << " " << final_state2[0]  << std::endl;

}
