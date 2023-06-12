#ifndef ILP_APPROXIMATORS_H
#define ILP_APPROXIMATORS_H

#include <map>
#include <tuple>
#include <vector>
#include <sstream>

#include <boost/random.hpp>
#include <boost/dynamic_bitset.hpp>

#include "allegro/definitions.h"
#include "allegro/decode_bitset.h"
#include "ortools/linear_solver/linear_solver.h"

std::vector<std::pair<std::string, std::string>> randomized_rounding(
    std::vector<operations_research::MPVariable *> &feasible_solutions,
    boost::dynamic_bitset<> all_containers_bitset,
    std::map<boost::dynamic_bitset<>, std::pair<char, boost::dynamic_bitset<>>> &coversets,
    std::size_t multiplicity,
    std::size_t num_trials,
    std::ostringstream &log_buffer);

#endif