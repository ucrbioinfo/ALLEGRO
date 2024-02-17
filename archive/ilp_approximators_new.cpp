#include <thread>
#include <fstream>
#include <stdlib.h>
#include <iostream>
#include <filesystem>

#include "allegro/definitions.h"
#include "allegro/decode_bitset.h"
#include "allegro/ilp_approximators.h"

#include "absl/time/time.h"
#include "ortools/base/logging.h"
#include "ortools/sat/cp_model.h"
#include "ortools/sat/cp_model.pb.h"
#include "ortools/sat/cp_model_solver.h"
#include "ortools/sat/model.h"
#include "ortools/sat/sat_parameters.pb.h"
#include "ortools/util/sorted_interval_list.h"

std::vector<GuideStruct> sat_solver(
    std::vector<operations_research::MPVariable *> &feasible_solutions,
    boost::dynamic_bitset<> all_containers_bitset,
    std::map<boost::dynamic_bitset<>, std::pair<double, boost::dynamic_bitset<>>> &coversets,
    std::size_t multiplicity,
    std::size_t beta,
    std::size_t early_stopping_patience_s,
    std::ostringstream &log_buffer
    )
{
    std::cout << BLUE << "> " << RESET << "Setting up the ILP problem..." << std::endl;
    // Recreate coversets using feasible variables
    std::map<boost::dynamic_bitset<>, std::pair<double, boost::dynamic_bitset<>>> all_feasible_coverings;

    for (auto var_ptr : feasible_solutions) {
        boost::dynamic_bitset<> bitset(var_ptr->name());

        double score = coversets[bitset].first;
        boost::dynamic_bitset<> species_hit_by_this_guide = coversets[bitset].second;

        all_feasible_coverings[bitset] = std::pair<double, boost::dynamic_bitset<>>(score, species_hit_by_this_guide);
    }

    // Create a model with the ILP CP-SAT backend.
    operations_research::sat::CpModelBuilder cp_model;

    // --------------------------------------------------
    // -------------- VARIABLE CREATION -----------------
    // --------------------------------------------------
    std::map<boost::dynamic_bitset<>, std::set<boost::dynamic_bitset<>>> hit_containers;

    // Maps a bitset (representing a protospace sequence) to an OR-TOOLS variable.
    std::map<boost::dynamic_bitset<>, operations_research::sat::BoolVar> map_seq_to_vars;
    std::vector<operations_research::sat::BoolVar> all_vars;

    auto it = all_feasible_coverings.begin();
    while (it != all_feasible_coverings.end())
    {
        double score = it->second.first;

        boost::dynamic_bitset<> guide_seq_bitset = it->first;
        boost::dynamic_bitset<> container_bitset = it->second.second;

        // Find the first container hit by this guide and while there are containers left to process...
        size_t set_bit_index = container_bitset.find_first();
        while (set_bit_index != boost::dynamic_bitset<>::npos)
        {
            // Create a onehot bitvector for this one species specifically.
            // For three species, this would be:
            //  first iteration: bit vector 001 for the first species,
            //  second iteration: 010 for the second, and then third iteration 100 for the third.
            boost::dynamic_bitset<> container_onehot(container_bitset.size());
            container_onehot.set(set_bit_index);

            // Indicate that this container is hit by this guide.
            // For example, 001: 00110011... when container 1 is hit by this ATAT... guide
            // and on the next iteration, if ATAT hits container 2 as well:
            // 010: 00110011...
            hit_containers[container_onehot].insert(guide_seq_bitset);

            // Find the next container hit by this guide for the next iteration.
            // Returns boost::dynamic_bitset<>::npos if no other bits are set.
            set_bit_index = container_bitset.find_next(set_bit_index);
        }

        std::string buffer;
        boost::to_string(guide_seq_bitset, buffer);  // Bitset to string: "0101010110111"...
        operations_research::sat::BoolVar const var = cp_model.NewBoolVar().WithName(buffer);
        map_seq_to_vars[guide_seq_bitset] = var;
        all_vars.push_back(var);

        it++;
    }

    // --------------------------------------------------
    // ------------- CONSTRAINT CREATION ----------------
    // --------------------------------------------------
    for (auto i : hit_containers)
    {
        std::vector<operations_research::sat::BoolVar> vars_for_this_container;
        std::set<boost::dynamic_bitset<>> seq_bitsets_set = i.second;

        for (auto j : seq_bitsets_set)
        {
            vars_for_this_container.push_back(map_seq_to_vars[j]);
        }

        cp_model.AddGreaterOrEqual(operations_research::sat::LinearExpr::Sum(vars_for_this_container), 1);
    }

    hit_containers.clear(); // Mark memory as free

    if (beta > 0)
    {
        cp_model.AddLessOrEqual(operations_research::sat::LinearExpr::Sum(all_vars), beta);
        cp_model.Maximize(operations_research::sat::LinearExpr::Sum(all_vars));
    }
    else {
        cp_model.Minimize(operations_research::sat::LinearExpr::Sum(all_vars));
    }

    // Solving part
    operations_research::sat::Model model;

    // Sets a time limit
    std::cout << BLUE << "> " << RESET << "Solving the ILP within the given time limit of " << early_stopping_patience_s << " seconds..." << std::endl;
    operations_research::sat::SatParameters parameters;
    parameters.set_max_time_in_seconds(early_stopping_patience_s);
    model.Add(NewSatParameters(parameters));

    // Solve
    const operations_research::sat::CpSolverResponse response = SolveCpModel(cp_model.Build(), &model);

    // Check that the problem has a solution
    log_buffer << "Status: " << response.status() << std::endl;
    if (response.status() == operations_research::sat::OPTIMAL)
    {
        // std::cout << BLUE << "> Status: " << response.status() << RESET << std::endl;
        std::cout << BLUE << "> " << RESET << "Found an " << BLUE << "optimal" << RESET << " solution!" << std::endl;
        log_buffer << "The ILP has an optimal solution!" << std::endl;
    }
    else if (response.status() == operations_research::sat::FEASIBLE)
    {
        // std::cout << BLUE << "> Status: " << response.status() << RESET << std::endl;
        std::cout << BLUE << "> " << RESET << "Found a " << BLUE << "feasible" << RESET << " solution. Increasing the patience may result in a better solution." << std::endl;
        log_buffer << "The ILP has a feasible solution." << std::endl;
    }
    else
    {
        std::cout << RED << "> Status: " << response.status() << RESET << std::endl;
        std::cout << RED << "> The ILP problem cannot be solved. Exiting." << RESET << std::endl;
        log_buffer << "The ILP problem cannot be solved." << std::endl;
        return std::vector<GuideStruct>();
    }

    std::vector<operations_research::sat::BoolVar> sat_feasible_solutions;

    for (auto i : map_seq_to_vars)
        {
            boost::dynamic_bitset<> seq_bitset = i.first;
            operations_research::sat::BoolVar var = i.second;

            if (operations_research::sat::SolutionIntegerValue(response, var) > 0.0)
            {
                log_buffer << decode_bitset(seq_bitset) << " with solution value: " << operations_research::sat::SolutionIntegerValue(response, var) << std::endl;

                sat_feasible_solutions.push_back(var);
            }
        }

        map_seq_to_vars.clear();

        std::size_t len_solutions = sat_feasible_solutions.size();
        std::cout << BLUE << "> " << RESET << "The final set consists of: " << len_solutions << " guides." << std::endl;

        std::vector<GuideStruct> decoded_winners;

        log_buffer << "The final set consists of:" << std::endl;
        for (auto var_ptr : sat_feasible_solutions)
        {
            boost::dynamic_bitset<> bitset(var_ptr.Name());

            double score = all_feasible_coverings[bitset].first;
            boost::dynamic_bitset<> species_hit_by_this_guide = all_feasible_coverings[bitset].second;

            std::string buffer;
            boost::to_string(species_hit_by_this_guide, buffer);

            std::string decoded_bitset = decode_bitset(var_ptr.Name());

            GuideStruct guide;
            guide.sequence = decoded_bitset;
            guide.score = score;
            guide.species_hit = buffer;

            decoded_winners.push_back(guide);

            log_buffer << decoded_bitset << std::endl;
        }

    log_buffer << "Size of the final set: " << decoded_winners.size() << std::endl;

    return decoded_winners;
}
