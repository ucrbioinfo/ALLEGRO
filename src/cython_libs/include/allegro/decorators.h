#ifndef DECORATORS_H
#define DECORATORS_H

#include <map>
#include <sstream>
#include <algorithm>

#include <boost/dynamic_bitset.hpp>

void decorate_with_monophonic(
    std::size_t cut_multiplicity,
    std::size_t monophonic_threshold,
    std::ostringstream &log_buffer,
    std::map<boost::dynamic_bitset<>, std::pair<double, boost::dynamic_bitset<>>> &coversets);

#endif