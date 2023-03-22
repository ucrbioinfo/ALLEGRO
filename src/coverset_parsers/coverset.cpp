#include <vector>
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
        this->all_species_bitset.set();
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
            boost::dynamic_bitset<> bitset(this->num_species);
            bitset.set(species_id);

            std::pair<unsigned char, boost::dynamic_bitset<>> p(score, bitset);

            this->coversets[encoded] = p;
        }
    }

    std::string CoversetCPP::decode_bitset(boost::dynamic_bitset<> &encoded)
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

    std::string CoversetCPP::decode_bitset(const std::string &encoded_str)
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
        std::cout << "Using randomized rounding with " << this->num_trials << " trials.\n";

        std::unordered_set<std::string> winners;
        std::size_t len_winners = this->num_species;
        std::size_t trial_with_smallest_size = 0;

        // On each invocation, dist(rng) returns a random floating-point value uniformly distributed in the range [min..max).
        boost::random::mt19937 rng(std::time(nullptr));
        boost::random::uniform_real_distribution<double> dist(0, 1);

        for (std::size_t trial = 1; trial < num_trials + 1; trial++)
        {
            // Species to cover
            boost::dynamic_bitset<> I_this_trial(this->num_species);
            std::unordered_set<std::string> winners_this_trial;
            // std::size_t iterations_this_trial = 0;

            while (!(I_this_trial & this->all_species_bitset).all())
            {
                // TODO This needs an exit condition. The loop can be infinite as of now.
                // iterations_this_trial++;
                for (auto var_ptr : feasible_solutions)
                {
                    operations_research::MPVariable *var = var_ptr;

                    if ((var->solution_value() == 1.0) || (var->solution_value() > dist(rng)))
                    {
                        boost::dynamic_bitset<> bitset(var->name());
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

        std::cout << "Winners are:\n";
        for (auto winner_str : winners)
        {
            boost::dynamic_bitset<> bitset(winner_str);
            boost::dynamic_bitset<> species_hit_by_this_guide = this->coversets[bitset].second;

            std::string buffer;
            boost::to_string(species_hit_by_this_guide, buffer);

            std::string decoded_bitset = decode_bitset(winner_str);
            decoded_winners.push_back(std::pair<std::string, std::string> (decoded_bitset, buffer));

            std::cout <<  decoded_bitset << std::endl;
        }

        return decoded_winners;
    }

    std::vector<std::pair<std::string, std::string>> CoversetCPP::ortools_solver()
    {
        // Create the linear solver with the GLOP backend.
        std::unique_ptr<operations_research::MPSolver> solver(operations_research::MPSolver::CreateSolver("GLOP"));

        std::unordered_map<boost::dynamic_bitset<>, std::unordered_set<boost::dynamic_bitset<>>> hit_species;

        std::unordered_set<boost::dynamic_bitset<>> species_already_hit_by_unique_guide;
        std::unordered_set<boost::dynamic_bitset<>> marked_for_death;

        for (auto i : this->coversets)
        {
            boost::dynamic_bitset<> guide_seq_bits = i.first;
            boost::dynamic_bitset<> species_bitset = i.second.second;

            // We want to keep only one guide per species where that guide hits
            //  only this species and none other.
            // If a species is already hit by a one-hitting-guide and we encounter
            //  another one, mark it for deletion
            //   from this->coversets and skip adding it for hit_species.
            if (species_bitset.count() == 1)
            {
                if ((species_already_hit_by_unique_guide.find(species_bitset) != species_already_hit_by_unique_guide.end())) {
                    marked_for_death.insert(guide_seq_bits);
                    continue;
                }
                else {
                    size_t set_bit_index = species_bitset.find_first();
                    boost::dynamic_bitset<> species_onehot(species_bitset.size());
                    species_onehot.set(set_bit_index);

                    hit_species[species_onehot].insert(guide_seq_bits);
                    species_already_hit_by_unique_guide.insert(species_bitset);

                    set_bit_index = species_bitset.find_next(set_bit_index);
                }
            }
            else
            {
                size_t set_bit_index = species_bitset.find_first();
                while (set_bit_index != boost::dynamic_bitset<>::npos)
                {
                    boost::dynamic_bitset<> species_onehot(species_bitset.size());
                    species_onehot.set(set_bit_index);

                    hit_species[species_onehot].insert(guide_seq_bits);

                    set_bit_index = species_bitset.find_next(set_bit_index);
                }
            }
        }

        // Remove redundant guides from further processing.
        // We do not want to make variables for these.
        for (auto i : marked_for_death)
        {
            this->coversets.erase(i);
        }

        marked_for_death.clear();
        species_already_hit_by_unique_guide.clear();

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

        hit_species.clear();

        LOG(INFO) << "Number of constraints = " << solver->NumConstraints();
        // --------------------------------------------------

        // Set the objective and solve.
        objective->SetMinimization();

        const operations_research::MPSolver::ResultStatus result_status = solver->Solve();

        // Check that the problem has an optimal solution.
        if (result_status != operations_research::MPSolver::OPTIMAL)
        {
            LOG(FATAL) << "The problem does not have an optimal solution!";
        }

        std::cout << "Status: " << result_status << std::endl;

        // Save the feasible variables.
        std::vector<operations_research::MPVariable *> feasible_solutions;
        for (auto i : map_seq_to_vars)
        {
            boost::dynamic_bitset<> seq_bitset = i.first;
            operations_research::MPVariable *var = i.second;

            if (var->solution_value() > 0.0)
            {
                std::cout << decode_bitset(seq_bitset) << " with solution value " << var->solution_value() << std::endl;
                feasible_solutions.push_back(var);
            }
        }

        map_seq_to_vars.clear();

        std::size_t len_solutions = feasible_solutions.size();
        std::cout << "Number of feasible candidate guides: " << len_solutions << std::endl;

        // --------------------------------------------------
        // -------------- RANDOMIZED ROUND ------------------
        // --------------------------------------------------
        if (len_solutions > 0) {
            return randomized_rounding(feasible_solutions);
        } else {
            return std::vector<std::pair<std::string, std::string>>();
        }
    }
}