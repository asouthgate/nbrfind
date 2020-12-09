#ifndef STATE_DATA_HPP
#define STATE_DATA_HPP

class StateData {
    public:
        StateData(int max_row_length);
        ~StateData() { delete[] state_arr; }
        void freeze(int prev_lower_bound, int prev_upper_bound, int d, int maxi);
        void swap_pointers();
        void init_state_array(int rowsize);
        void init_state_quintuple(int len1, int len2);
        void print_arrays();
        void print_debug();
        
        int MAX_ROW_SIZE;
        // TODO: make private with public interface
        // Fundamental arrays used for calculations:
        // L is aln dist; M is n matches
        // NM is number of mismatches I believe
        // NN is number of number of pairs with an N
        int* L0; int* L1; int* L2;
        int* M0; int* M1; int* M2;
        int* NM0; int* NM1; int* NM2;
        int* NN0; int* NN1; int* NN2;

        int h;
        int lower_bound;
        int upper_bound;
        int d_resume = 0;
        int maxi_resume = 0;

        int* state_arr;
        int state_quintuple[5] = {0,0,0,0,0};

    private:
};

#endif
