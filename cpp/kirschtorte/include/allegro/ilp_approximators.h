#ifndef ILP_APPROXIMATORS_H
#define ILP_APPROXIMATORS_H

#include <map>
#include <tuple>
#include <vector>
#include <sstream>

#include <boost/dynamic_bitset.hpp>

#include "allegro/definitions.h"
#include "allegro/guide_struct.h"
#include "allegro/decode_bitset.h"
#include "ortools/linear_solver/linear_solver.h"

void sat_solver(
    std::vector<operations_research::MPVariable *> &feasible_solutions,
    std::map<boost::dynamic_bitset<>, std::pair<double, boost::dynamic_bitset<>>> &coversets,
    std::size_t const &multiplicity,
    std::size_t &beta,
    std::size_t &early_stopping_patience,
    bool const &enable_solver_diagnostics,
    std::string const &output_directory,
    std::ostringstream &log_buffer,
    std::vector<GuideStruct> &solution_set);

#endif