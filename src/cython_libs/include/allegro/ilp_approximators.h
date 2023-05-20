#ifndef ILP_APPROXIMATORS_H
#define ILP_APPROXIMATORS_H

#include <vector>
#include <sstream>
#include <unordered_map>
#include <unordered_set>

#include <boost/random.hpp>
#include <boost/dynamic_bitset.hpp>

#include "allegro/decode_bitset.h"
#include "ortools/linear_solver/linear_solver.h"

std::vector<std::pair<std::string, std::string>> randomized_rounding(
    std::vector<operations_research::MPVariable *> &feasible_solutions,
    boost::dynamic_bitset<> &all_species_bitset,
    std::unordered_map<boost::dynamic_bitset<>, std::pair<char, boost::dynamic_bitset<>>> &coversets,
    std::size_t num_containers,
    std::size_t num_trials,
    std::ostringstream &log_buffer);

#endif