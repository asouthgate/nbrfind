#include <vector>
#include <string>

const std::string five_prime = "ATGGAGAGCCTTGTCCCTGG";
const std::string three_prime = "TTAACTTTAATCTCACATAG";

bool trim_to_cds(std::string& s);

void convert_nonstandard_to_N(std::string& s);

bool filter(const std::string& s, int min_length, int max_length);

std::vector<std::pair<std::string, std::string>> read_fasta(std::string fasta_fname, int min_len, int max_len);

