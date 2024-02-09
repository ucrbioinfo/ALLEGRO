#ifndef ORTOOLSSOLVER_H
#define ORTOOLSSOLVER_H

#include <string>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>

#include <boost/dynamic_bitset.hpp>
#include "ortools/linear_solver/linear_solver.h"

namespace solvers
{
    class ORToolsSolver
    {
    private:
        std::size_t num_trials;
        std::size_t num_species;
        std::size_t guide_length;
        boost::dynamic_bitset<> all_species_bitset;

        std::string decode_bitset(const std::string &encoded_str);
        std::string decode_bitset(boost::dynamic_bitset<> &encoded);
        std::unordered_map<boost::dynamic_bitset<>, std::pair<unsigned char, boost::dynamic_bitset<>>> *coversets;

        std::unordered_set<std::string> randomized_rounding(
            std::vector<operations_research::MPVariable *> feasible_solutions);

    public:
        ORToolsSolver(
            std::size_t num_species,
            std::size_t guide_length,
            std::size_t num_trials,
            std::unordered_map<boost::dynamic_bitset<>, std::pair<unsigned char, boost::dynamic_bitset<>>> &coversets);

        ~ORToolsSolver();
        std::unordered_set<std::string> ortools_solver();
    };
}

#endif