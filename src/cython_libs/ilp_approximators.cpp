#include <thread>
#include <fstream>
#include <iostream>
#include <filesystem>

#include "allegro/definitions.h"
#include "allegro/decode_bitset.h"
#include "allegro/ilp_approximators.h"

#include "absl/time/time.h"
#include "ortools/linear_solver/linear_solver.h"

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

    // Create a linear solver with the ILP SAT backend.
    std::unique_ptr<operations_research::MPSolver> solver(operations_research::MPSolver::CreateSolver("SAT"));
    
    // solver->EnableOutput();
    
    int time_limit_ms = early_stopping_patience_s * 1000; // Time limit in milliseconds
    absl::Duration time_limit = absl::Milliseconds(time_limit_ms); // Convert milliseconds to absl::Duration
    solver->SetTimeLimit(time_limit);

    const double infinity = solver->infinity(); // Used for constraints to denote >= 1 and <= 1

    // May return 0 when not able to detect
    const auto processor_count = std::thread::hardware_concurrency();

    if (processor_count > 1) {
        absl::Status status = solver->SetNumThreads(processor_count);

        if (!status.ok()) {
            std::cout << RED << "> Error setting number of threads: " << RESET << status.message() << std::endl;
        }
    }

    // --------------------------------------------------
    // -------------- VARIABLE CREATION -----------------
    // --------------------------------------------------
    operations_research::MPObjective *const objective = solver->MutableObjective();
    operations_research::MPConstraint *beta_constraint;

    if (beta > 0)
    {
        beta_constraint = solver->MakeRowConstraint(-infinity, beta);
    }

    std::map<boost::dynamic_bitset<>, std::set<boost::dynamic_bitset<>>> hit_containers;

    // Maps a bitset (representing a protospace sequence) to an OR-TOOLS variable.
    std::map<boost::dynamic_bitset<>, operations_research::MPVariable *> map_seq_to_vars;

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
        boost::to_string(guide_seq_bitset, buffer);
        operations_research::MPVariable *const var = solver->MakeIntVar(0, 1, buffer);

        map_seq_to_vars[guide_seq_bitset] = var;

        if (beta > 0)
        {
            objective->SetCoefficient(var, score);
            beta_constraint->SetCoefficient(var, 1);
        }
        else
        {
            objective->SetCoefficient(var, 1);
        }

        it++;
    }

    // --------------------------------------------------
    // ------------- CONSTRAINT CREATION ----------------
    // --------------------------------------------------
    for (auto i : hit_containers)
    {
        std::vector<operations_research::MPVariable *> vars_for_this_container;
        std::set<boost::dynamic_bitset<>> seq_bitsets_set = i.second;

        for (auto j : seq_bitsets_set)
        {
            vars_for_this_container.push_back(map_seq_to_vars[j]);
        }

        // All species/genes must be covered by at least multiplicity guide(s)
        operations_research::MPConstraint *const constraint = solver->MakeRowConstraint(multiplicity, infinity);
        for (auto k : vars_for_this_container)
        {
            constraint->SetCoefficient(k, 1);
        }
    }

    hit_containers.clear(); // Mark memory as free
    
    // Set the appropriate objective.
    if (beta > 0)
    {
        // LOG(INFO) << "Beta: " << beta << " - Maximizing the the set size considering beta." << std::endl;
        // this->log_buffer << "Beta: " << beta << " - Maximizing the the set size considering beta." << std::endl;

        objective->SetMaximization();
    }
    else
    {
        // LOG(INFO) << "Beta: " << beta << " - Minimizing the the set size. Disregarding scores and beta." << std::endl;
        // this->log_buffer << "Beta: " << beta << " - Minimizing the the set size. Disregarding scores and beta." << std::endl;

        objective->SetMinimization();
    }

    // Solve the integer linear program
    std::cout << BLUE << "> " << RESET << "Solving the ILP within the given time limit of " << early_stopping_patience_s << " seconds..." << std::endl;
    const operations_research::MPSolver::ResultStatus result_status = solver->Solve();

    // Check that the problem has a solution.
    log_buffer << "Status: " << result_status << std::endl;
    if (result_status == operations_research::MPSolver::OPTIMAL)
    {
        std::cout << BLUE << "> Status: " << result_status << RESET << std::endl;
        std::cout << BLUE << "> " << RESET << "The ILP has an optimal solution!" << std::endl;
        log_buffer << "The ILP has an optimal solution!" << std::endl;
    }
    else if (result_status == operations_research::MPSolver::FEASIBLE)
    {
        std::cout << BLUE << "> Status: " << result_status << RESET << std::endl;
        std::cout << BLUE << "> " << RESET << "The ILP has a feasible solution." << std::endl;
        log_buffer << "The ILP has a feasible solution." << std::endl;
    }
    else
    {
        std::cout << RED << "> Status: " << result_status << RESET << std::endl;
        std::cout << RED << "> The ILP cannot be solved. Exiting." << RESET << std::endl;
        log_buffer << "The ILP cannot be solved." << std::endl;
        return std::vector<GuideStruct>();
    }

    std::vector<operations_research::MPVariable *> sat_feasible_solutions;

    for (auto i : map_seq_to_vars)
        {
            boost::dynamic_bitset<> seq_bitset = i.first;
            operations_research::MPVariable *var = i.second;

            if (var->solution_value() > 0.0)
            {
                log_buffer << decode_bitset(seq_bitset) << " with solution value: " << var->solution_value() << std::endl;

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
            boost::dynamic_bitset<> bitset(var_ptr->name());

            double score = all_feasible_coverings[bitset].first;
            boost::dynamic_bitset<> species_hit_by_this_guide = all_feasible_coverings[bitset].second;

            std::string buffer;
            boost::to_string(species_hit_by_this_guide, buffer);

            std::string decoded_bitset = decode_bitset(var_ptr->name());

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


// FUNCTION IS DEPRECATED - ARCHIVAL CODE
// FUNCTION IS DEPRECATED - ARCHIVAL CODE
// FUNCTION IS DEPRECATED - ARCHIVAL CODE

// std::vector<GuideStruct> randomized_rounding(
//     std::vector<operations_research::MPVariable *> &feasible_solutions,
//     boost::dynamic_bitset<> all_containers_bitset,
//     std::map<boost::dynamic_bitset<>, std::pair<double, boost::dynamic_bitset<>>> &coversets,
//     std::size_t multiplicity,
//     std::size_t num_trials,
//     std::ostringstream &log_buffer,
//     std::string output_directory
//     )
// {

//     std::filesystem::path dirpath(output_directory);
//     std::filesystem::path filepath = dirpath / "trial_sizes.csv";

//     std::ofstream log_file(filepath);

//     log_file << "trial,library_size" << std::endl;

//     std::cout << BLUE << "> " << RESET << "Using randomized rounding with " << num_trials << " trials." << std::endl;
//     log_buffer << "Using randomized rounding with " << num_trials << " trials." << std::endl;

//     std::size_t original_size = all_containers_bitset.size();
//     std::size_t new_size = original_size * multiplicity;

//     boost::dynamic_bitset<> extended_all_containers_bitset(new_size);
//     all_containers_bitset.resize(new_size);

//     for (std::size_t i = 0; i < multiplicity; i++)
//     {
//         extended_all_containers_bitset <<= original_size;
//         extended_all_containers_bitset |= all_containers_bitset;
//     }

//     std::set<std::string> winners;
//     std::size_t len_winners = INT_MAX; // Infinity, initially
//     std::size_t trial_with_smallest_size = 0;

//     // On each invocation, dist(rng) returns a random floating-point value
//     // uniformly distributed in the range [min, max).
//     boost::random::mt19937 rng(std::time(nullptr));
//     boost::random::uniform_real_distribution<double> dist(0, 1);

//     std::size_t trial = 1;
//     while (trial < num_trials + 1)
//     {
//         if (trial % 1000 == 0)
//         {
//             // Writes on the same line
//             std::cout << BLUE << "\r> " << RESET << "Trial " << trial << std::flush;
//         }

//         // Guide containers to cover
//         boost::dynamic_bitset<> I_this_trial(new_size);
//         std::set<std::string> winners_this_trial;
//         std::size_t iterations_this_trial = 100000; // while-loop exit condition in case of bad luck

//         while ((I_this_trial != extended_all_containers_bitset) && (iterations_this_trial != 0))
//         {
//             iterations_this_trial--;

//             for (auto var_ptr : feasible_solutions)
//             {
//                 operations_research::MPVariable *var = var_ptr;

//                 if (winners_this_trial.find(var->name()) != winners_this_trial.end())
//                 {
//                     continue;
//                 }

//                 if ((var->solution_value() == 1.0) || (var->solution_value() > dist(rng)))
//                 {
//                     // Encoded binary DNA sequence
//                     boost::dynamic_bitset<> bitset(var->name());

//                     // The species bit vector hit by this binary DNA sequence
//                     boost::dynamic_bitset<> species_hit_by_this_guide = coversets[bitset].second;

//                     species_hit_by_this_guide.resize(new_size);
//                     boost::dynamic_bitset<> old_I(new_size);

//                     while (species_hit_by_this_guide.find_first() != boost::dynamic_bitset<>::npos)
//                     {
//                         old_I = I_this_trial;
//                         I_this_trial |= species_hit_by_this_guide;
//                         species_hit_by_this_guide &= old_I;
//                         species_hit_by_this_guide <<= original_size;
//                     }

//                     winners_this_trial.insert(var->name());
//                 }
//             }
//         }

//         // double score_sum_this_trial = 0.0;
//         // for (auto itr = winners_this_trial.begin(); itr != winners_this_trial.end(); itr++)
//         // {
//         //     boost::dynamic_bitset<> bitset(*itr);
//         //     score_sum_this_trial += this->coversets[bitset].first;
//         // }

//         std::size_t len_winners_this_trial = winners_this_trial.size();


//         // Log the size of the library for this trial
//         log_file << trial << "," << len_winners_this_trial << std::endl;

//         if (len_winners_this_trial <= len_winners)
//         {
//             winners = winners_this_trial;
//             len_winners = len_winners_this_trial;
//             trial_with_smallest_size = trial;
//         }

//         trial += 1;
//     }

//     if (trial > 1000)
//     {
//         std::cout << std::endl;
//     }

//     std::vector<GuideStruct> decoded_winners;

//     // std::cout << BLUE << "> " << RESET << "Winners are:" << std::endl;
//     log_buffer << "Winners are:" << std::endl;
//     for (auto winner_str : winners)
//     {
//         boost::dynamic_bitset<> bitset(winner_str);

//         double score = coversets[bitset].first;
//         boost::dynamic_bitset<> species_hit_by_this_guide = coversets[bitset].second;

//         std::string buffer;
//         boost::to_string(species_hit_by_this_guide, buffer);

//         std::string decoded_bitset = decode_bitset(winner_str);

//         GuideStruct guide;
//         guide.sequence = decoded_bitset;
//         guide.score = score;
//         guide.species_hit = buffer;

//         decoded_winners.push_back(guide);

//         // std::cout << BLUE << "> " << RESET << decoded_bitset << std::endl;
//         log_buffer << decoded_bitset << std::endl;
//     }

//     log_file.close();

//     return decoded_winners;
// }