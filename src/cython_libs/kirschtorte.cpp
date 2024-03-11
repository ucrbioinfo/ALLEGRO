#include <iostream>

#include "allegro/logger.h"
#include "allegro/decorators.h"
#include "allegro/kirschtorte.h"
#include "allegro/definitions.h"
#include "allegro/decode_bitset.h"
#include "allegro/ilp_approximators.h"

#include "absl/time/time.h"
#include "ortools/linear_solver/linear_solver.h"

namespace Kirschtorte
{
    Kirschtorte::Kirschtorte(
        std::size_t num_containers,
        std::size_t guide_length,
        std::size_t early_stopping_patience_s,
        std::string output_directory,
        bool enable_solver_diagnostics)
    {
        this->num_containers = num_containers; // A container may be a gene, species, or chromosome.
        this->guide_length = guide_length;     // Twenty (20) for cas9.
        this->early_stopping_patience_s = early_stopping_patience_s;
        this->output_directory = output_directory;
        this->bits_required_to_store_seq = guide_length * 2; // Each nucleotide A/C/T/G in a guide
                                                             // can be represented by 2 bits.
        this->enable_solver_diagnostics = enable_solver_diagnostics;
    }

    Kirschtorte::~Kirschtorte() {}

    int Kirschtorte::encode_and_save_dna(
        std::string seq,
        double score,
        std::size_t container_id)
    {
        std::string encoded_str("");

        for (auto base : seq)
        {
            switch (base)
            {
            case 'A':
                encoded_str += sA_SHIFT;
                break;
            case 'C':
                encoded_str += sC_SHIFT;
                break;
            case 'G':
                encoded_str += sG_SHIFT;
                break;
            case 'T':
                encoded_str += sT_SHIFT;
                break;
            default:
                std::cerr << "Invalid DNA character: '" << base << "' encountered in " << seq << " in container id " << container_id << ". Skipping sequence.\n";
                return 1;
            }
        }

        // This is the encoded bitset for a guide's sequence
        boost::dynamic_bitset<> encoded_bitset = boost::dynamic_bitset<>(encoded_str);

        // If the sequence already exists in the memory
        if (this->coversets.find(encoded_bitset) != this->coversets.end())
        {
            // Set the appropriate bit to indicate that this new container is hit by this guide.
            this->coversets[encoded_bitset].second.set(container_id);

            // TODO: Update the average score of this guide
            // need to keep track of how many times this guide has been seen: n
            // new_average = old_average + (new_item - old_average) / (n + 1)
        }
        else
        {
            // If it is the first time we encounter this sequence
            // Create a bitset with width equal to the number of containers
            // For example, if we have 10 species, we will have 0000000000
            boost::dynamic_bitset<> bitset(this->num_containers);
            bitset.set(container_id); // Set the appropriate bit for this container to 1

            // Indicate that this guide (its encoded bitset) has a pair<score, and the container it targets>
            this->coversets[encoded_bitset] = std::pair<double, boost::dynamic_bitset<>>(score, bitset);
        }

        return 0;
    }

    std::vector<GuideStruct> Kirschtorte::setup_and_solve(
        std::size_t monophonic_threshold,
        std::size_t multiplicity,
        std::size_t beta,
        std::size_t seed_length,
        std::size_t mismatched_allowed_after_seed,
        bool precluster)
    {
        beta = 0;
        
        ::google::InitGoogleLogging("ALLEGRO VON AMIR");

        // Create a linear solver with the GLOP backend.
        std::unique_ptr<operations_research::MPSolver> solver(operations_research::MPSolver::CreateSolver("GLOP"));
        const double infinity = solver->infinity(); // Used for constraints to denote >= 1 and <= 1

        // --------------------------------------------------
        // ----------------- DECORATORS ---------------------
        // --------------------------------------------------
        if (monophonic_threshold > 0)
        {
            // Within this function, guides deemed not needed have their scores set to 0.
            // These guides will be removed below where we have: if (score <= 0) {...}
            decorate_with_monophonic(multiplicity, monophonic_threshold, this->log_buffer, this->coversets);
        }

        if (precluster)
        {
            if (multiplicity == 1)
            {
                decorate_with_clustering(
                    seed_length,
                    multiplicity,
                    mismatched_allowed_after_seed,
                    this->log_buffer,
                    this->coversets);
            }
            else
            {
                decorate_with_clustering_multiplicity(
                    seed_length,
                    multiplicity,
                    mismatched_allowed_after_seed,
                    this->log_buffer,
                    this->coversets);
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

        // Maps a bitset (representing a target container (gene or species) to
        // a set of bitsets (representing a target sequence).
        std::map<boost::dynamic_bitset<>, std::set<boost::dynamic_bitset<>>> hit_containers;

        // Maps a bitset (representing a target sequence) to an OR-TOOLS variable.
        std::map<boost::dynamic_bitset<>, operations_research::MPVariable *> map_seq_to_vars;

        auto it = this->coversets.begin();
        while (it != this->coversets.end())
        {
            double score = it->second.first;

            // Remove useless guides.
            // Either predicted inefficient by the scorer, or marked as not needed above.
            if (score <= 0)
            {
                it = this->coversets.erase(it);
                continue;
            }

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
            operations_research::MPVariable *const var = solver->MakeNumVar(0.0, 1.0, buffer);

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

        LOG(INFO) << "Number of variables (guides): " << solver->NumVariables();
        this->log_buffer << "Number of variables (guides): " << solver->NumVariables() << std::endl;
        // --------------------------------------------------

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

        if (beta > 0)
        {
            LOG(INFO) << "Number of constraints + 1 (beta): " << solver->NumConstraints();
            this->log_buffer << "Number of constraints + 1 (beta): " << solver->NumConstraints() << std::endl;
        }
        else
        {
            LOG(INFO) << "Number of constraints: " << solver->NumConstraints();
            this->log_buffer << "Number of constraints: " << solver->NumConstraints() << std::endl;
        }
        // --------------------------------------------------

        LOG(INFO) << "Cut multiplicity: " << multiplicity << std::endl;
        this->log_buffer << "Cut multiplicity: " << multiplicity << std::endl;

        // Set the appropriate objective.
        if (beta > 0)
        {
            LOG(INFO) << "Beta: " << beta << " - Maximizing the the set size considering beta." << std::endl;
            this->log_buffer << "Beta: " << beta << " - Maximizing the the set size considering beta." << std::endl;

            objective->SetMaximization();
        }
        else
        {
            LOG(INFO) << "Beta: " << beta << " - Minimizing the the set size. Disregarding scores and beta." << std::endl;
            this->log_buffer << "Beta: " << beta << " - Minimizing the the set size. Disregarding scores and beta." << std::endl;

            objective->SetMinimization();
        }

        // Solve the linear program
        operations_research::MPSolver::ResultStatus result_status = solver->Solve();

        // Check that the problem has a solution.
        if (result_status == operations_research::MPSolver::OPTIMAL)
        {
            // std::cout << BLUE << "> Status: " << result_status << RESET << std::endl;
            std::cout << BLUE << "> " << RESET << "The LP problem has an " << BLUE << "optimal" << RESET << " solution!" << std::endl;
            this->log_buffer << "The LP problem has an optimal solution." << std::endl;
        }
        else if (result_status == operations_research::MPSolver::FEASIBLE)
        {
            // std::cout << BLUE << "> Status: " << result_status << RESET << std::endl;
            std::cout << BLUE << "> " << RESET << "The LP problem has a " << BLUE << "feasible" << RESET << " solution." << std::endl;
            this->log_buffer << "The LP problem has a feasible solution." << std::endl;
        }
        else
        {
            std::cout << RED << "> Status: " << result_status << RESET << std::endl;
            std::cout << RED << "> The LP problem cannot be solved." << RESET << std::endl;
            this->log_buffer << "The LP problem cannot be solved. Status: " << result_status << std::endl;

            if (this->enable_solver_diagnostics)
            {
                std::size_t counter = 1;
                bool fixed_beta = false;

                std::cout << BLUE << "> " << RESET << "Diagnosing constraints by iteratively relaxing them and resolving..." << std::endl;
                
                for (auto constraint : solver->constraints())
                {
                    // Temporarily relax the constraint
                    std::cout << BLUE "\r> " << RESET << "Relaxing constraint " << counter << "/" << solver->NumConstraints() << "..." << std::flush;
                    constraint->SetBounds(-infinity, infinity);

                    result_status = solver->Solve();  // Resolve with relaxed constraint.
                                                    // If the constraint is BETA, this is easy and fast to solve.

                    // Check feasibility
                    if ((result_status == operations_research::MPSolver::OPTIMAL) || (result_status == operations_research::MPSolver::FEASIBLE))
                    {
                        std::cout << BLUE "\n> " << RESET << "Relaxing constraint " << constraint->name() << " makes the problem feasible." << std::endl;
                        this->log_buffer << "Relaxing constraint " << constraint->name() << " makes the problem feasible." << std::endl;

                        // Was Beta the bad constraint? Binary search for the best Beta
                        if (constraint->name() == "BETA")
                        {
                            std::cout << BLUE "> " << RESET << "Looking for the lowest feasible beta..." << std::endl;
                            
                            // Remove all scores from guides and resolve with minimization.
                            // The objective value of the minimization will be the new (smallest) beta.
                            for (auto var : solver->variables())
                            {
                                objective->SetCoefficient(var, 1);
                            }
                            
                            objective->SetMinimization();

                            result_status = solver->Solve();  // Resolve with relaxed constraint

                            if ((result_status == operations_research::MPSolver::OPTIMAL) || (result_status == operations_research::MPSolver::FEASIBLE))
                            {
                                std::size_t min_beta = objective->Value() + 1;
                                // The actual objective value is a FLOAT value.
                                // truncating the fractional part results in infeasible. One workaround is to add 1 to fix it.
                                
                                std::cout << BLUE "> " << RESET << "Increasing beta to " << min_beta << " makes the problem feasible. This may increase if we need to solve the ILP." << std::endl;
                                log_buffer << "Increasing beta " << min_beta << " makes the problem feasible. This may increase if we need to solve the ILP." << std::endl;
                                
                                beta = min_beta;
                                fixed_beta = true;
                                break;
                            }
                        }
                        else
                        {
                            std::cout << BLUE "> " << RESET << "Exiting." << std::endl; 
                            log_info(this->log_buffer, this->output_directory);
                            return std::vector<GuideStruct>();
                        }
                    }
                    counter++;
                }

                if (fixed_beta == false)
                {
                    std::cout << RED << "> Unfortunately ALLEGRO could not find the issue. Possibly more than a single constrait is defective. The LP problem cannot be solved. You can try iteratively removing genes and/or species and resolving." << RESET << std::endl;
                    this->log_buffer << "Relaxing each constraint and resolving did not find the issue. The LP problem cannot be solved. Exiting." << std::endl;
                    log_info(this->log_buffer, this->output_directory);
                    return std::vector<GuideStruct>();
                }
            }
            else
            {
                std::cout << RED << "> Exiting. Enable diagnostics in config.yaml to iteratively look for a possibly bad constraint." << RESET << std::endl;
                log_info(this->log_buffer, this->output_directory);
                return std::vector<GuideStruct>();
            }
        }

        // Save the feasible variables.
        // A feasible variable is any variable with a solution value greater than 0.
        // It has a chance to be in the final solution.
        std::vector<operations_research::MPVariable *> feasible_solutions;
        std::size_t num_have_fractional_vars = 0;

        for (auto i : map_seq_to_vars)
        {
            boost::dynamic_bitset<> seq_bitset = i.first;
            operations_research::MPVariable *var = i.second;

            if (var->solution_value() > 0.0)
            {
                if (var->solution_value() < 1.0)
                {
                    num_have_fractional_vars++;
                }

                this->log_buffer << decode_bitset(seq_bitset) << " with solution value: " << var->solution_value() << std::endl;

                feasible_solutions.push_back(var);
            }
        }

        map_seq_to_vars.clear();

        std::size_t len_solutions = feasible_solutions.size();
        // --------------------------------------------------
        // ---------------------- ILP -----------------------
        // --------------------------------------------------

        // vector of guide structs [guide sequence, score, container bit vector that it cuts]
        std::vector<GuideStruct> solution_set;

        if (len_solutions > 0)
        {
            if (num_have_fractional_vars > 0)
            {
                std::cout << BLUE << "> " << RESET << "Number of feasible candidate guides: " << len_solutions << "." << std::endl;
                std::cout << BLUE << "> " << RESET << "Switching to the ILP solver as " << num_have_fractional_vars << " residual guides remain." << std::endl;
                this->log_buffer << "Number of feasible candidate guides: " << len_solutions << "." << std::endl;
                
                // SAT solver with time limit
                solution_set = sat_solver(
                    feasible_solutions,
                    this->coversets,
                    multiplicity,
                    beta,
                    this->early_stopping_patience_s,
                    this->enable_solver_diagnostics,
                    this->output_directory,
                    this->log_buffer);
            }
            else
            {
                std::cout << BLUE << "> " << RESET <<  "The final set consists of " << len_solutions << " guides." << std::endl;
                this->log_buffer << "The final set consists of:" << std::endl;

                for (auto var_ptr : feasible_solutions)
                {
                    boost::dynamic_bitset<> bitset(var_ptr->name());

                    double score = this->coversets[bitset].first;
                    boost::dynamic_bitset<> species_hit_by_this_guide = this->coversets[bitset].second;

                    std::string buffer;
                    boost::to_string(species_hit_by_this_guide, buffer);

                    std::string decoded_bitset = decode_bitset(var_ptr->name());

                    GuideStruct guide;
                    guide.sequence = decoded_bitset;  // This will be in binary and NOT "ACTGTG..."
                    guide.score = score;
                    guide.species_hit = buffer;

                    solution_set.push_back(guide);

                    log_buffer << decoded_bitset << std::endl;
                }
            }
        }
        else
        {
            std::cout << "No feasible solutions." << len_solutions << std::endl;
            this->log_buffer << "No feasible solutions." << len_solutions << std::endl;
        }

        // Write all buffer messages to disk
        log_info(this->log_buffer, this->output_directory);

        return solution_set;
    }
}