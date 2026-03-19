#ifndef DECORATORS_H
#define DECORATORS_H

#include <map>
#include <set>
#include <sstream>
#include <algorithm>
#include <boost/dynamic_bitset.hpp>

void decorate_with_monophonic(
    std::size_t cut_multiplicity,
    std::size_t monophonic_threshold,
    std::ostringstream &log_buffer,
    std::map<boost::dynamic_bitset<>, std::pair<double, boost::dynamic_bitset<>>> &coversets);

void decorate_with_clustering(
    std::size_t seed_length,
    std::size_t cut_multiplicity,
    std::size_t mismatched_allowed_after_seed,
    std::ostringstream &log_buffer,
    std::map<boost::dynamic_bitset<>, std::pair<double, boost::dynamic_bitset<>>> &coversets);

void decorate_with_clustering_multiplicity(
    std::size_t seed_length,
    std::size_t cut_multiplicity,
    std::size_t mismatched_allowed_after_seed,
    std::ostringstream &log_buffer,
    std::map<boost::dynamic_bitset<>, std::pair<double, boost::dynamic_bitset<>>> &coversets);

#endif