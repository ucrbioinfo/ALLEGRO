#include <vector>
#include <iostream>
#include <unordered_set>

#include "allegro/coverset.h"
#include "ortools/linear_solver/linear_solver.h"


#define bA_SHIFT 0b00
#define bC_SHIFT 0b01
#define bG_SHIFT 0b10
#define bT_SHIFT 0b11

#define sA_SHIFT "00"
#define sC_SHIFT "01"
#define sG_SHIFT "10"
#define sT_SHIFT "11"


namespace coversets {
    CoversetCPP::CoversetCPP(std::size_t num_species, std::size_t guide_length) {
        this->num_species = num_species;
        this->guide_length = guide_length;
    }


    CoversetCPP::~CoversetCPP() {}


    void CoversetCPP::encode_and_save_dna(
        std::string &seq,
        unsigned char score,
        unsigned short species_id)
    {
        std::size_t bits_required_to_store_seq = this->guide_length * 2;

        boost::dynamic_bitset<> encoded(bits_required_to_store_seq);

        boost::dynamic_bitset<> A_shift(bits_required_to_store_seq, bA_SHIFT);
        boost::dynamic_bitset<> C_shift(bits_required_to_store_seq, bC_SHIFT);
        boost::dynamic_bitset<> G_shift(bits_required_to_store_seq, bG_SHIFT);
        boost::dynamic_bitset<> T_shift(bits_required_to_store_seq, bT_SHIFT);

        for (size_t i = 0; i < this->guide_length; i++) {
            // Shift the bitset to the left by 2 bits
            encoded <<= 2;

            // Set the last two bits based on the current DNA character
            switch (seq[i]) {
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
        if (coversets.find(encoded) != coversets.end()) {
            this->coversets[encoded].second[species_id] = 1;
        }
        else {
            boost::dynamic_bitset<> bitset(this->num_species);

            bitset.set(species_id);

            std::pair<unsigned char, boost::dynamic_bitset<>> pair(score, bitset);

            this->coversets[encoded] = pair;
        }
    }

    // For DEBUGGING -- Decodes and returns the first bitset in the coversets object
    std::string CoversetCPP::get_str() {
        boost::dynamic_bitset<> encoded = this->coversets.begin()->first;
        std::string decoded = "";

        for (std::size_t i = 0; i < this->guide_length; i++) {

            // Take the 2 least sigbits
            boost::dynamic_bitset<> LSBs(this->guide_length * 2, 0b11);
            LSBs &= encoded;

            std::string buffer;
            boost::to_string(LSBs, buffer);  // Insert into buffer as a string

            std::string last_two_chars = buffer.substr(buffer.length() - 2);

            // Check what the last two bits were
            if (last_two_chars == sA_SHIFT) {
                decoded.push_back('A');
            }
            else if (last_two_chars == sC_SHIFT) {
                decoded.push_back('C');
            }
            else if (last_two_chars == sG_SHIFT) {
                decoded.push_back('G');
            }
            else if (last_two_chars == sT_SHIFT) {
                decoded.push_back('T');
            }

            encoded >>= 2;
        }

        std::reverse(decoded.begin(), decoded.end());
        std::cout << decoded << std::endl;
        return decoded;
    }

    void CoversetCPP::ortools_solver() {
        // // Create the linear solver with the GLOP backend.
        // std::unique_ptr<operations_research::MPSolver> solver(operations_research::MPSolver::CreateSolver("GLOP"));

        // // std::unordered_map<std::bitset<40>,  std::pair<unsigned char, std::bitset<339>>> coversets;

        // std::unordered_map<std::bitset<339>, std::unordered_set<std::bitset<40>>> hit_species;

        // for (auto it = this->coversets.begin(); it != this->coversets.end(); it++) {
        //     std::bitset<40> key = it->first;
        //     std::pair<unsigned char, std::bitset<339>> value = it->second;

        //     unsigned char score = value.first;
        //     std::bitset<339> species_bitset = value.second;



        // }


        // const double infinity = solver->infinity();
        // operations_research::MPObjective* const objective = solver->MutableObjective();        
        // std::unordered_map<std::bitset<40>, operations_research::MPVariable *> vars;

        // for (auto it = this->coversets.begin(); it != this->coversets.end(); it++) {
        //     std::bitset<40> seq = it->first;

        //     operations_research::MPVariable* const var = solver->MakeNumVar(0.0, 1, seq.to_string());
            
        //     vars[seq] = var;
        //     objective->SetCoefficient(var, 1);
        // }

        // LOG(INFO) << "Number of variables = " << solver->NumVariables();


        // for (auto it = hit_species.begin(); it != hit_species.end(); it++) {
        //     std::vector<operations_research::MPVariable*> vars_for_this_species;
        //     std::unordered_set<std::bitset<40>> seqs = it->second;

        //     for (auto it = seqs.begin(); it != seqs.end(); it++) {
        //         vars_for_this_species.push_back(vars[*it]);
        //     }
            
        //     operations_research::MPConstraint* const constraint = solver->MakeRowConstraint(0.0, infinity);

        //     for (auto it = vars_for_this_species.begin(); it != vars_for_this_species.end(); it++) {
        //         constraint->SetCoefficient(*it, 1);
        //     }
        // }

        // LOG(INFO) << "Number of constraints = " << solver->NumConstraints();
        
        // // // Free up memory
        // // this->coversets.clear();

        // objective->SetMinimization();

        // std::string lol = "idk";
        // solver->ExportModelAsLpFormat(false, &lol);

        // const operations_research::MPSolver::ResultStatus result_status = solver->Solve();

        // // Check that the problem has an optimal solution.
        // if (result_status != operations_research::MPSolver::OPTIMAL) {
        // LOG(FATAL) << "The problem does not have an optimal solution!";
        // }

        // for (auto it = vars.begin(); it != vars.end(); it++) {

        // }

        // LOG(INFO) << "Solution:" << std::endl;
        // LOG(INFO) << "Objective value = " << objective->Value();
    }
}