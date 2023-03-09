#ifndef COVERSETCPP_H
#define COVERSETCPP_H

#include <string>
#include <algorithm>
#include <unordered_map>
#include <boost/dynamic_bitset.hpp>


namespace coversets {
    class CoversetCPP {
        public:
            std::size_t num_species;
            std::size_t guide_length;
            std::unordered_map<boost::dynamic_bitset<>, std::pair<unsigned char, boost::dynamic_bitset<>>> coversets;
            
            CoversetCPP(std::size_t num_species, std::size_t guide_length);
            ~CoversetCPP();

            void ortools_solver();
            std::string get_str();
            void encode_and_save_dna(std::string& seq, unsigned char score, unsigned short species_id);
    };
}

#endif