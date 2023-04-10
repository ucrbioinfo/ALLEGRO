#include <vector>
#include <fstream>
#include <iostream>
#include <boost/random.hpp>

#include "allegro/coverset.h"

#define bA_SHIFT 0b00
#define bC_SHIFT 0b01
#define bG_SHIFT 0b10
#define bT_SHIFT 0b11

#define sA_SHIFT "00"
#define sC_SHIFT "01"
#define sG_SHIFT "10"
#define sT_SHIFT "11"

std::string generate_log_filename()
{
    int count = 1;
    std::string filename;

    while (true)
    {
        // Generate a filename based on the current count
        if (count == 1)
        {
            filename = "last_runs_log.txt";
        }
        else
        {
            filename = "last_runs_log_" + std::to_string(count) + ".txt";
        }

        // Check if the file exists
        std::ifstream file(filename);
        if (!file)
        {
            return filename;
        }

        count++;
    }
}

void log_info(std::ostringstream& log_buffer)
{
    std::ofstream log_file(generate_log_filename());

    if (!log_file.is_open())
    {
        std::cerr << "Unable to open log file!" << std::endl;
    }

    log_file << log_buffer.str();
    log_file.close();
}

namespace coversets
{
    CoversetCPP::CoversetCPP(
        std::size_t num_species,
        std::size_t guide_length,
        std::size_t num_trials)
    {
        this->num_species = num_species;
        this->guide_length = guide_length;
        this->num_trials = num_trials;
        this->bits_required_to_store_seq = guide_length * 2;

        this->all_species_bitset = boost::dynamic_bitset<>(num_species);
    }

    CoversetCPP::~CoversetCPP() {}

    void CoversetCPP::encode_and_save_dna(
        std::string &seq,
        unsigned char score,
        unsigned short species_id)
    {
        // Assuming we are using alphabet of 4 A/C/T/G, each is represented by 2 bits.
        boost::dynamic_bitset<> encoded(this->bits_required_to_store_seq);

        boost::dynamic_bitset<> A_shift(this->bits_required_to_store_seq, bA_SHIFT);
        boost::dynamic_bitset<> C_shift(this->bits_required_to_store_seq, bC_SHIFT);
        boost::dynamic_bitset<> G_shift(this->bits_required_to_store_seq, bG_SHIFT);
        boost::dynamic_bitset<> T_shift(this->bits_required_to_store_seq, bT_SHIFT);

        for (size_t i = 0; i < guide_length; i++)
        {
            // Shift the bitset to the left by 2 bits
            encoded <<= 2;

            // Set the last two bits based on the current DNA character
            switch (seq[i])
            {
            case 'A':
                encoded |= A_shift;
                break;
            case 'C':
                encoded |= C_shift;
                break;
            case 'G':
                encoded |= G_shift;
                break;
            case 'T':
                encoded |= T_shift;
                break;
            default:
                std::cerr << "Invalid DNA character: " << seq[i] << '\n';
            }
        }

        // If seq already exists
        if (coversets.find(encoded) != coversets.end())
        {
            // Set the appropriate bit to indicate which species is hit by this guide.
            this->coversets[encoded].second.set(species_id);
        }
        else
        {
            // If it is the first time we encounter this sequence
            boost::dynamic_bitset<> bitset(this->num_species);
            bitset.set(species_id);

            this->coversets[encoded] = std::pair<unsigned char, boost::dynamic_bitset<>>(score, bitset);
        }

        // Keep a record of which species should be hit. We compare against this later in randomized_rounding
        //  to determine if a set of feasible solutions hits all the required species or not.
        // We do this record keeping here instead of in the constructor because some species
        //  may not contain any guides and should be excluded from consideration (so we don't set those bits).
        this->all_species_bitset.set(species_id);
    }

    std::string CoversetCPP::decode_bitset(boost::dynamic_bitset<> encoded)
    {
        std::string decoded = "";
        std::string buffer;
        std::string last_two_chars;

        for (std::size_t i = 0; i < this->guide_length; i++)
        {
            // Take the 2 least sigbits
            boost::dynamic_bitset<> LSBs(encoded.size(), 0b11);
            LSBs &= encoded;

            boost::to_string(LSBs, buffer); // Insert into buffer as a string

            last_two_chars = buffer.substr(buffer.length() - 2);

            // Check what the last two bits were
            if (last_two_chars == sA_SHIFT)
            {
                decoded.push_back('A');
            }
            else if (last_two_chars == sC_SHIFT)
            {
                decoded.push_back('C');
            }
            else if (last_two_chars == sG_SHIFT)
            {
                decoded.push_back('G');
            }
            else if (last_two_chars == sT_SHIFT)
            {
                decoded.push_back('T');
            }

            encoded >>= 2;
        }

        std::reverse(decoded.begin(), decoded.end());
        return decoded;
    }

    std::string CoversetCPP::decode_bitset(const std::string encoded_str)
    {
        boost::dynamic_bitset<> encoded(encoded_str);

        std::string decoded = "";
        std::string buffer;
        std::string last_two_chars;

        for (std::size_t i = 0; i < this->guide_length; i++)
        {
            // Take the 2 least sigbits
            boost::dynamic_bitset<> LSBs(encoded.size(), 0b11);
            LSBs &= encoded;

            boost::to_string(LSBs, buffer); // Insert into buffer as a string

            last_two_chars = buffer.substr(buffer.length() - 2);

            // Check what the last two bits were
            if (last_two_chars == sA_SHIFT)
            {
                decoded.push_back('A');
            }
            else if (last_two_chars == sC_SHIFT)
            {
                decoded.push_back('C');
            }
            else if (last_two_chars == sG_SHIFT)
            {
                decoded.push_back('G');
            }
            else if (last_two_chars == sT_SHIFT)
            {
                decoded.push_back('T');
            }

            encoded >>= 2;
        }

        std::reverse(decoded.begin(), decoded.end());
        return decoded;
    }

    // For DEBUGGING -- Decodes and returns the first bitset in the coversets object
    std::string CoversetCPP::get_str()
    {
        boost::dynamic_bitset<> encoded = this->coversets.begin()->first;
        std::string decoded = "";

        for (std::size_t i = 0; i < this->guide_length; i++)
        {
            // Take the 2 least sigbits
            boost::dynamic_bitset<> LSBs(this->bits_required_to_store_seq, 0b11);
            LSBs &= encoded;

            std::string buffer;
            boost::to_string(LSBs, buffer); // Insert into buffer as a string

            std::string last_two_chars = buffer.substr(buffer.length() - 2);

            // Check what the last two bits were
            if (last_two_chars == sA_SHIFT)
            {
                decoded.push_back('A');
            }
            else if (last_two_chars == sC_SHIFT)
            {
                decoded.push_back('C');
            }
            else if (last_two_chars == sG_SHIFT)
            {
                decoded.push_back('G');
            }
            else if (last_two_chars == sT_SHIFT)
            {
                decoded.push_back('T');
            }

            encoded >>= 2;
        }

        std::reverse(decoded.begin(), decoded.end());
        std::cout << decoded << std::endl;
        return decoded;
    }

    std::vector<std::pair<std::string, std::string>> CoversetCPP::randomized_rounding(
        std::vector<operations_research::MPVariable *> feasible_solutions)
    {
        // Algorithm source:
        // https://web.archive.org/web/20230325165759/https://theory.stanford.edu/~trevisan/cs261/lecture08.pdf
        std::cout << "Using randomized rounding with " << this->num_trials << " trials." << std::endl;
        this->log_buffer << "Using randomized rounding with " << this->num_trials << " trials." << std::endl;

        std::unordered_set<std::string> winners;
        std::size_t len_winners = this->num_species;
        std::size_t trial_with_smallest_size = 0;

        // On each invocation, dist(rng) returns a random floating-point value
        //  uniformly distributed in the range [min, max).
        boost::random::mt19937 rng(std::time(nullptr));
        boost::random::uniform_real_distribution<double> dist(0, 1);

        for (std::size_t trial = 1; trial < num_trials + 1; trial++)
        {
            // Species to cover
            boost::dynamic_bitset<> I_this_trial(this->num_species);
            std::unordered_set<std::string> winners_this_trial;
            std::size_t iterations_this_trial = 0; // while-loop exit condition in case of bad luck

            while ((I_this_trial != this->all_species_bitset) && (iterations_this_trial != 100000))
            {
                iterations_this_trial++;

                for (auto var_ptr : feasible_solutions)
                {
                    operations_research::MPVariable *var = var_ptr;

                    if ((var->solution_value() == 1.0) || (var->solution_value() > dist(rng)))
                    {
                        // Encoded binary DNA sequence
                        boost::dynamic_bitset<> bitset(var->name());

                        // The species bit vector hit by this binary DNA sequence
                        boost::dynamic_bitset<> species_hit_by_this_guide = this->coversets[bitset].second;

                        I_this_trial |= species_hit_by_this_guide;

                        winners_this_trial.insert(var->name());
                    }
                }
            }

            // double score_sum_this_trial = 0.0;
            // for (auto itr = winners_this_trial.begin(); itr != winners_this_trial.end(); itr++)
            // {
            //     boost::dynamic_bitset<> bitset(*itr);
            //     score_sum_this_trial += this->coversets[bitset].first;
            // }

            std::size_t len_winners_this_trial = winners_this_trial.size();

            if (len_winners_this_trial <= len_winners)
            {
                winners = winners_this_trial;
                len_winners = len_winners_this_trial;
                trial_with_smallest_size = trial;
            }
        }

        std::vector<std::pair<std::string, std::string>> decoded_winners;

        std::cout << "Winners are:" << std::endl;
        this->log_buffer << "Winners are:" << std::endl;
        for (auto winner_str : winners)
        {
            boost::dynamic_bitset<> bitset(winner_str);
            boost::dynamic_bitset<> species_hit_by_this_guide = this->coversets[bitset].second;

            std::string buffer;
            boost::to_string(species_hit_by_this_guide, buffer);

            std::string decoded_bitset = decode_bitset(winner_str);
            decoded_winners.push_back(std::pair<std::string, std::string>(decoded_bitset, buffer));

            std::cout << decoded_bitset << std::endl;
            this->log_buffer << decoded_bitset << std::endl;
        }

        return decoded_winners;
    }

    std::vector<std::pair<std::string, std::string>> CoversetCPP::ortools_solver()
    {
        // Create the linear solver with the GLOP backend.
        std::unique_ptr<operations_research::MPSolver> solver(operations_research::MPSolver::CreateSolver("GLOP"));

        std::unordered_map<boost::dynamic_bitset<>, std::unordered_set<boost::dynamic_bitset<>>> hit_species;

        std::unordered_set<boost::dynamic_bitset<>> species_already_hit_by_unique_guide;
        std::unordered_set<boost::dynamic_bitset<>> zerstoeren;

        for (auto i : this->coversets)
        {
            boost::dynamic_bitset<> guide_seq_bits = i.first;

            unsigned char score = i.second.first;
            boost::dynamic_bitset<> species_bitset = i.second.second;

            // Below, we want to keep only one guide per species where that guide hits
            //  only this species and none other.
            // If a species is already hit by a one-hitting-guide and we encounter
            //  another one, mark it for deletion from this->coversets and
            //   skip adding it for hit_species.
            if (species_bitset.count() <= 3)
            {
                // and if this species already has a representative guide that hits it...
                if ((species_already_hit_by_unique_guide.find(species_bitset) != species_already_hit_by_unique_guide.end()))
                {
                    // The new guide is not needed.
                    // TODO May compare scores here and below and replace if better.
                    // Mark the new guide for deletion and carry on.
                    zerstoeren.insert(guide_seq_bits);
                    continue;
                }
                // If this species still needs a representative guide...
                else
                {
                    // Indicate that this species has a representative guide now and does not need another one later.
                    species_already_hit_by_unique_guide.insert(species_bitset);
                }
            }
            // Find the first species hit by this guide and while there are species left to process...
            size_t set_bit_index = species_bitset.find_first();
            while (set_bit_index != boost::dynamic_bitset<>::npos)
            {
                // Create a onehot bitvector for this one species specifically.
                // For three species, this would be:
                //  first iteration: bit vector 001 for the first species,
                //  second iteration: 010 for the second, and then third iteration 100 for the third.
                boost::dynamic_bitset<> species_onehot(species_bitset.size());
                species_onehot.set(set_bit_index);

                // Indicate that this species is hit by this guide.
                // For example, 001: 00110011... when species 1 is hit by this ATAT... guide
                // and on the next iteration, if ATAT hits species 2 as well:
                //  010: 00110011...
                hit_species[species_onehot].insert(guide_seq_bits);

                // Find the next species hit by this guide for the next iteration.
                // Returns boost::dynamic_bitset<>::npos if no other bits are set.
                set_bit_index = species_bitset.find_next(set_bit_index);
            }
        }

        // Space saving: Remove redundant guides from further processing.
        // We do not want to make solver variables for these.
        for (auto i : zerstoeren)
        {
            this->coversets.erase(i);
        }

        zerstoeren.clear();                          // Mark memory as free
        species_already_hit_by_unique_guide.clear(); // Mark memory as free

        // --------------------------------------------------
        // -------------- VARIABLE CREATION -----------------
        // --------------------------------------------------
        operations_research::MPObjective *const objective = solver->MutableObjective();
        std::unordered_map<boost::dynamic_bitset<>, operations_research::MPVariable *> map_seq_to_vars;

        for (auto i : this->coversets)
        {
            boost::dynamic_bitset<> seq_bitset = i.first;
            unsigned char score = i.second.first;

            std::string buffer;
            boost::to_string(seq_bitset, buffer);
            operations_research::MPVariable *const var = solver->MakeNumVar(0.0, 1, buffer);

            map_seq_to_vars[seq_bitset] = var;
            objective->SetCoefficient(var, score);
        }

        LOG(INFO) << "Number of variables = " << solver->NumVariables();
        this->log_buffer << "Number of variables = " << solver->NumVariables() << std::endl;
        // --------------------------------------------------

        // --------------------------------------------------
        // ------------- CONSTRAINT CREATION ----------------
        // --------------------------------------------------
        const double infinity = solver->infinity();
        for (auto i : hit_species)
        {
            std::vector<operations_research::MPVariable *> vars_for_this_species;
            std::unordered_set<boost::dynamic_bitset<>> seq_bitsets_set = i.second;

            for (auto j : seq_bitsets_set)
            {
                vars_for_this_species.push_back(map_seq_to_vars[j]);
            }

            operations_research::MPConstraint *const constraint = solver->MakeRowConstraint(1.0, infinity);
            for (auto k : vars_for_this_species)
            {
                constraint->SetCoefficient(k, 1);
            }
        }

        hit_species.clear(); // Mark memory as free

        LOG(INFO) << "Number of constraints (species) = " << solver->NumConstraints();
        this->log_buffer << "Number of constraints (species) = " << solver->NumConstraints() << std::endl;
        // --------------------------------------------------

        // Set the objective and solve.
        objective->SetMinimization();

        const operations_research::MPSolver::ResultStatus result_status = solver->Solve();

        // Check that the problem has an optimal solution.
        std::cout << "Status: " << result_status << std::endl;
        this->log_buffer << "Status: " << result_status << std::endl;
        if (result_status != operations_research::MPSolver::OPTIMAL)
        {
            LOG(FATAL) << "The problem does not have an optimal solution!";
            this->log_buffer << "The problem does not have an optimal solution!" << std::endl;
        }

        // Save the feasible variables.
        // A feasible variable is any variable with a solution value greater than 0.
        // It has a chance to be in the final solution.
        std::vector<operations_research::MPVariable *> feasible_solutions;
        for (auto i : map_seq_to_vars)
        {
            boost::dynamic_bitset<> seq_bitset = i.first;
            operations_research::MPVariable *var = i.second;

            if (var->solution_value() > 0.0)
            {
                std::cout << decode_bitset(seq_bitset) << " with solution value " << var->solution_value() << std::endl;
                this->log_buffer << decode_bitset(seq_bitset) << " with solution value " << var->solution_value() << std::endl;

                feasible_solutions.push_back(var);
            }
        }

        map_seq_to_vars.clear();

        std::size_t len_solutions = feasible_solutions.size();
        std::cout << "Number of feasible candidate guides: " << len_solutions << std::endl;
        this->log_buffer << "Number of feasible candidate guides: " << len_solutions << std::endl;
        // --------------------------------------------------
        // -------------- RANDOMIZED ROUND ------------------
        // --------------------------------------------------
        log_info(this->log_buffer);

        if (len_solutions > 0)
        {
            return randomized_rounding(feasible_solutions);
        }
        else
        {
            // Empty -- A problem with 0 feasible solutions (empty inputs or no guides in the fasta files)
            //  still returns an OPTIMAL status by GLOP. Return an empty vector in this edge case.
            // Why would you input no guides? :/
            return std::vector<std::pair<std::string, std::string>>();
        }
    }
}