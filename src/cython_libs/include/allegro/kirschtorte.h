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
            std::size_t num_containers,
            std::size_t guide_length,
            std::size_t early_stopping_patience_s,
            std::string output_directory,
            bool enable_solver_diagnostics);

        ~Kirschtorte();

        std::vector<GuideStruct> setup_and_solve(
            std::size_t monophonic_threshold,
            std::size_t multiplicity,
            std::size_t beta,
            std::size_t seed_length,
            std::size_t mismatched_allowed_after_seed,
            bool precluster);

        int encode_and_save_dna(
            std::string seq,
            double score,
            std::size_t container_id);

    private:
        std::size_t early_stopping_patience_s;
        std::size_t num_containers;
        std::size_t guide_length;
        std::string output_directory;
        std::size_t bits_required_to_store_seq;
        boost::dynamic_bitset<> all_containers_bitset;
        std::ostringstream log_buffer;
        bool enable_solver_diagnostics;

        // guide bitvector --> (score, bitvector of species it hits)
        // The width of the guide bitvector is 2x guide length AKA this->(std::size_t) bits_required_to_store_seq
        std::map<boost::dynamic_bitset<>, std::pair<double, boost::dynamic_bitset<>>> coversets;

        void create_lp_constraints(
            std::unique_ptr<operations_research::MPSolver> &solver,
            std::map<boost::dynamic_bitset<>, std::set<boost::dynamic_bitset<>>> const &hit_containers,
            std::map<boost::dynamic_bitset<>, operations_research::MPVariable *> &map_seq_to_vars,
            std::size_t const &multiplicity);

        void get_feasible_solutions(
            std::map<boost::dynamic_bitset<>, operations_research::MPVariable *> const &map_seq_to_vars,
            std::vector<operations_research::MPVariable *> &feasible_solutions,
            std::size_t &num_have_fractional_vars);

        operations_research::MPSolver::ResultStatus setup_and_solve_with_beta(
            std::unique_ptr<operations_research::MPSolver> const &solver,
            operations_research::MPObjective *const &objective,
            std::size_t const &beta,
            std::map<operations_research::MPVariable *, std::size_t> const &map_var_to_score);

        operations_research::MPSolver::ResultStatus setup_and_solve_without_beta(
            std::unique_ptr<operations_research::MPSolver> const &solver,
            operations_research::MPObjective *const &objective);

        std::string diagnose_lp(
            operations_research::MPSolver::ResultStatus result_status,
            bool const &enable_solver_diagnostics,
            std::unique_ptr<operations_research::MPSolver> const &solver,
            operations_research::MPObjective *const &objective,
            std::size_t &beta,
            std::ostringstream &log_buffer,
            std::string const &output_directory);
    };
}

#endif