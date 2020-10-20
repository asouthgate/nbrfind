class StateData {
    public:
        StateData(int* input_state_quintuple, int* input_state_arr);
        void freeze(int prev_lower_bound, int prev_upper_bound, int d, int maxi);
        void swap_pointers();
        
        // TODO: make private with public interface
        // Fundamental arrays used for calculations:
        // L is aln dist; M is n matches
        // NM is number of mismatches I believe
        // NN is number of number of pairs with an N
        int* L0; int* L1; int* L2;
        int* M0; int* M1; int* M2;
        int* NM0; int* NM1; int* NM2;
        int* NN0; int* NN1; int* NN2;
        int& h;
        int& lower_bound;
        int& upper_bound;
        int& d_resume;
        int& maxi_resume;

    private:
        int* state_arr;
        int* state_quintuple;
}
