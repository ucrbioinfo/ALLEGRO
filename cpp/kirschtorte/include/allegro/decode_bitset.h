#ifndef DECODE_BITSET_H
#define DECODE_BITSET_H

#include <boost/dynamic_bitset.hpp>

std::string decode_bitset(boost::dynamic_bitset<> encoded);
std::string decode_bitset(const std::string &encoded_str);

#endif