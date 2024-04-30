#include <thread>
#include <fstream>
#include <iostream>
#include <filesystem>

#include "allegro/logger.h"
#include "allegro/definitions.h"
#include "allegro/decode_bitset.h"
#include "allegro/ilp_approximators.h"

#include "absl/time/time.h"
#include "ortools/linear_solver/linear_solver.h"

void sat_solver(
    std::vector<operations_research::MPVariable *> &feasible_solutions,
    std::map<boost::dynamic_bitset<>, std::pair<double, boost::dynamic_bitset<>>> &coversets,
    std::size_t const &multiplicity,
    std::size_t &beta,
    std::size_t &early_stopping_patience,
    bool const &enable_solver_diagnostics,
    std::string const &output_directory,
    std::ostringstream &log_buffer,
    std::vector<GuideStruct> &solution_set)
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
    
    int time_limit = early_stopping_patience;
    absl::Duration time_limit_ms = absl::Milliseconds(time_limit * 1000);  // Convert milliseconds to absl::Duration
    solver->SetTimeLimit(time_limit_ms);

    const double infinity = solver->infinity();  // Used for constraints to denote >= 1 and <= 1

    // May return 0 when not able to detect
    const auto processor_count = std::thread::hardware_concurrency();
    if (processor_count > 1) {
        absl::Status status = solver->SetNumThreads(processor_count - 1);

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
        beta_constraint = solver->MakeRowConstraint(-infinity, beta, "BETA");
    }

    std::map<boost::dynamic_bitset<>, std::set<boost::dynamic_bitset<>>> hit_containers;

    // Maps a bitset (representing a target sequence) to an OR-TOOLS variable.
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
        objective->SetMaximization();
    }
    else
    {
        objective->SetMinimization();
    }

    // Solve the integer linear program
    std::cout << BLUE << "> " << RESET << "Solving the ILP within the given time limit of " << early_stopping_patience << " seconds..." << std::endl;
    operations_research::MPSolver::ResultStatus result_status = solver->Solve();
    
    bool solved = false;
    std::size_t num_retries = 1;
    std::size_t auto_time_lim = 120;

    while (!solved && (num_retries < 10))
    {
        if (result_status == operations_research::MPSolver::OPTIMAL)
        {
            std::cout << BLUE << "> " << RESET << "Found an " << BLUE << "optimal" << RESET << " solution!" << std::endl;
            log_buffer << "The ILP has an optimal solution!" << std::endl;
            solved = true;
            break;
        }
        else if (result_status == operations_research::MPSolver::FEASIBLE)
        {
            std::cout << BLUE << "> " << RESET << "Found a " << BLUE << "feasible" << RESET << " solution. Increasing the patience may result in finding an optimal solution." << std::endl;
            log_buffer << "The ILP has a feasible solution." << std::endl;
            solved = true;
            break;
        }
        // Time limit issue - NOT_SOLVED means not **yet** solved
        else if (result_status == operations_research::MPSolver::NOT_SOLVED)
        {
            std::cout << RED << "> The ILP problem cannot be solved within the given time limit of " << time_limit << " seconds." << RESET << std::endl;

            if (enable_solver_diagnostics)
            {
                if (early_stopping_patience < auto_time_lim)
                {
                    time_limit = auto_time_lim;
                }
                else
                {
                    time_limit = early_stopping_patience * 1.5;
                }

                time_limit_ms = absl::Milliseconds(time_limit * 1000);
                solver->SetTimeLimit(time_limit_ms);
                early_stopping_patience = time_limit;
                
                std::cout << BLUE << "> " << RESET << "Retrying with a time limit of " << early_stopping_patience << " seconds (attempt " << num_retries << "/10)..." << std::endl;

                if (beta > 0)
                {
                    float increment_percent = 0.01 * num_retries;
                    beta = beta * increment_percent + beta;
                    beta_constraint->SetBounds(-infinity, beta);
                    std::cout << BLUE << "> " << RESET << "Increasing beta by a factor of " << increment_percent << " to " << beta << "." << std::endl;
                }

                result_status = solver->Solve();
                // Go back to the while(!solved..) loop and check result_status again

                if ((result_status != operations_research::MPSolver::FEASIBLE) && (result_status != operations_research::MPSolver::OPTIMAL))
                {
                    num_retries++;
                }
            }
            else
            {
                std::cout << RED << "> Exiting. Retrying with a larger patience may result in a solved problem." << RESET << std::endl;
                log_info(log_buffer, output_directory);
                return;
            }
        }
        else
        {
            std::cout << RED << "> Status: " << result_status << RESET << std::endl;
            std::cout << RED << "> The ILP problem cannot be solved." << RESET << std::endl;
            log_buffer << "The ILP problem cannot be solved." << std::endl;

            if (enable_solver_diagnostics)
            {
                std::size_t counter = 1;
                bool fixed_beta = false;

                std::cout << BLUE << "> " << RESET << "Diagnosing constraints by iteratively relaxing them and resolving..." << std::endl;
                
                for (auto constraint : solver->constraints())
                    {
                        // Temporarily relax the constraint
                        std::cout << BLUE "\r> " << RESET << "Relaxing constraint " << counter << "/" << solver->NumConstraints() << "..." << std::flush;
                        constraint->SetBounds(-infinity, infinity);

                        time_limit_ms = absl::Milliseconds(60 * 1000);  // Convert milliseconds to absl::Duration
                        solver->SetTimeLimit(time_limit_ms);

                        result_status = solver->Solve();  // Resolve with relaxed constraint

                        // Check feasibility
                        if ((result_status == operations_research::MPSolver::OPTIMAL) || (result_status == operations_research::MPSolver::FEASIBLE))
                        {
                            std::cout << BLUE "\n> " << RESET << "Relaxing constraint " << constraint->name() << " makes the problem feasible." << std::endl;
                            log_buffer << "Relaxing constraint " << constraint->name() << " makes the problem feasible." << std::endl;

                            // Was Beta the bad constraint? Search for the best Beta
                            if (constraint->name() == "BETA")
                            {
                                std::cout << BLUE "> " << RESET << "Looking for lowest feasible beta..." << std::endl;

                                // Remove all scores from guides and resolve with minimization.
                                // The objective value of the minimization will be the new (smallest) beta.
                                for (auto var : solver->variables())
                                {
                                    objective->SetCoefficient(var, 1);
                                }
                                
                                objective->SetMinimization();

                                time_limit = early_stopping_patience;  // Time limit in milliseconds
                                time_limit_ms = absl::Milliseconds(time_limit * 1000);  // Convert milliseconds to absl::Duration
                                solver->SetTimeLimit(time_limit_ms);

                                result_status = solver->Solve();  // Resolve with relaxed constraint

                                if ((result_status == operations_research::MPSolver::OPTIMAL) || (result_status == operations_research::MPSolver::FEASIBLE))
                                {
                                    std::size_t min_beta = objective->Value();

                                    std::cout << BLUE "> " << RESET << "Increasing beta to " << min_beta << " makes the problem feasible." << std::endl;
                                    log_buffer << "Increasing beta " << min_beta << " makes the problem feasible." << std::endl;
                                    
                                    beta = min_beta;
                                    fixed_beta = true;
                                    break;
                                }
                            }
                            else
                            {
                                std::cout << BLUE "> " << RESET << "Exiting." << std::endl; 
                                log_info(log_buffer, output_directory);
                                return;
                            }
                        }
                        counter++;
                    }

                if (fixed_beta == false)
                {
                    std::cout << RED << "> Unfortunately ALLEGRO could not find the issue. Possibly more than a single constraint is defective. The ILP problem cannot be solved. You can try iteratively removing genes and/or species and resolving." << RESET << std::endl;
                    log_buffer << "Relaxing each constraint and resolving did not find the issue. The ILP problem cannot be solved. Exiting." << std::endl;
                    log_info(log_buffer, output_directory);
                    return;
                }
            }
            else
            {
                std::cout << RED << "> Exiting. Enable diagnostics in config.yaml to iteratively look for a possibly bad constraint." << RESET << std::endl;
                log_info(log_buffer, output_directory);
                return;
            }
        }
    }

    if (!solved)
    {
        std::cout << RED << "> Exiting after 10 failed diagnostic retries." << RESET << std::endl;
        log_info(log_buffer, output_directory);
        return;
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

        solution_set.push_back(guide);

        log_buffer << decoded_bitset << std::endl;
    }

    log_buffer << "Size of the final set: " << solution_set.size() << std::endl;
}
