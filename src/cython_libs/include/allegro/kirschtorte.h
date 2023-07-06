#ifndef KIRSCHTORTE_H
#define KIRSCHTORTE_H

#include <map>
#include <string>
#include <vector>
#include <sstream>
#include <boost/dynamic_bitset.hpp>

#include "allegro/guide_struct.h"

namespace Kirschtorte
{
    class Kirschtorte
    {
    public:
        Kirschtorte(
            std::size_t num_containers,
            std::size_t guide_length,
            std::size_t num_trials,
            std::string output_directory);

        ~Kirschtorte();

        std::vector<GuideStruct> setup_and_solve(
            std::size_t monophonic_threshold,
            std::size_t cut_multiplicity,
            std::size_t beta);

        int encode_and_save_dna(
            std::string seq,
            double score,
            std::size_t container_id);

    private:
        std::size_t num_trials;
        std::size_t num_containers;
        std::size_t guide_length;
        std::string output_directory;
        std::size_t bits_required_to_store_seq;
        boost::dynamic_bitset<> all_containers_bitset;
        std::ostringstream log_buffer;

        // guide bitvector --> (score, bitvector of species it hits)
        // The width of the guide bitvector is 2 * guide length AKA this->(std::size_t) bits_required_to_store_seq
        std::map<boost::dynamic_bitset<>, std::pair<double, boost::dynamic_bitset<>>> coversets;
    };
}

#endif