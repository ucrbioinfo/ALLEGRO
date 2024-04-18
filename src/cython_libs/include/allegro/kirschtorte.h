#ifndef KIRSCHTORTE_H
#define KIRSCHTORTE_H

#include <map>
#include <string>
#include <vector>
#include <sstream>
#include <algorithm>
#include <boost/dynamic_bitset.hpp>
#include "allegro/guide_struct.h"
#include "ortools/linear_solver/linear_solver.h"

namespace Kirschtorte
{
    class Kirschtorte
    {
    public:
        Kirschtorte(
            bool precluster,
            std::size_t beta,
            std::size_t seed_length,
            std::size_t guide_length,
            std::size_t multiplicity,
            std::size_t num_containers,
            std::size_t monophonic_threshold,
            std::size_t early_stopping_patience,
            std::size_t mismatched_allowed_after_seed,
            bool enable_solver_diagnostics,
            std::string output_directory);

        ~Kirschtorte();

        std::vector<GuideStruct> setup_and_solve();

        int encode_and_save_dna(
            std::string seq,
            double score,
            std::size_t container_id);

    private:
        bool precluster;
        std::size_t beta;
        std::size_t seed_length;
        std::size_t guide_length;
        std::size_t multiplicity;
        std::size_t num_containers;
        std::size_t monophonic_threshold;
        std::size_t early_stopping_patience;
        std::size_t mismatched_allowed_after_seed;
        bool enable_solver_diagnostics;
        std::string output_directory;

        // Used to keep a rolling average of guides' scores
        double old_score;
        std::map<boost::dynamic_bitset<>, std::size_t> num_times_guide_was_seen;

        double infinity;
        operations_research::MPSolver *solver;
        operations_research::MPObjective *objective;

        std::ostringstream log_buffer;
        std::size_t bits_required_to_store_seq;
        boost::dynamic_bitset<> all_containers_bitset;
        
        // Maps an ortools var object to its predicted score, or 1.
        std::map<operations_research::MPVariable *, std::size_t> map_var_to_score;

        // Maps a bitset (representing a target container (gene or species) to a set of bitsets (representing a target sequence).
        std::map<boost::dynamic_bitset<>, std::set<boost::dynamic_bitset<>>> hit_containers;

        // Maps a bitset (representing a target sequence) to an OR-TOOLS variable.
        std::map<boost::dynamic_bitset<>, operations_research::MPVariable *> map_seq_to_vars;


        // guide bitvector --> (score, bitvector of species it hits)
        // The width of the guide bitvector is 2x guide length AKA this->(std::size_t) bits_required_to_store_seq
        std::map<boost::dynamic_bitset<>, std::pair<double, boost::dynamic_bitset<>>> coversets;

        void create_lp_constraints();
        operations_research::MPSolver::ResultStatus setup_and_solve_without_beta();
        operations_research::MPSolver::ResultStatus setup_and_solve_with_beta(std::size_t const &beta);

        std::string diagnose_lp(
            operations_research::MPSolver::ResultStatus result_status,
            std::size_t &beta);

        std::string fix_beta(std::size_t &beta);

        void get_feasible_solutions(
            std::vector<operations_research::MPVariable *> &feasible_solutions,
            std::size_t &num_have_fractional_vars);

    };
}

#endif