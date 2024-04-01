#include "allegro/decorators.h"
// MP THRESHOLD ALGORITHM.
// This function modifies the &coverset variable. Its primary purpose is to reduce the linear search space.
// Thus reducing the memory requirement by reducing the required number of variables in the LP.
//
// It directly removes redundant guides where it can, and adjusts the scores of other redundant guides
// to 0 when they cannot be directly removed.
//
// Redundant guides with a score of 0 will be removed in kirschtorte.cpp where coversets is iterated over again.
// The idea is that combinatorially selected container subsets of length mp_threshold need only have cut_multiplicity
// number of guides to represent them.
//
// For example, if you have 3 species where each species has 10 guides, mp_threshold 1, 
// cut_multiplicity 1, then each species need have only cut_multiplicity=1 best mp(1)-hitter guide*
// to represent it. We may discard the other 9 * 3 = 27 guides to save memory in the LP formulation phase.
// 
// *An mp(1)-hitter guide is a guide that hits ONLY a single species and NONE other simultaneously.
// We may safely discard this single-hitting guide without tampering with the linear program solution.
// If there exists a guide (with any score) that covers all 3 species at once, it will not be affected
// by this algorithm since in this example, it is not an mp(1)-hitting guide, it is a 3-hitting guide.
void decorate_with_monophonic(
    std::size_t cut_multiplicity,  // aka m guides required per container
    std::size_t monophonic_threshold, // aka combinatorially selected container subsets of length mp_threshold
    std::ostringstream &log_buffer,
    std::map<boost::dynamic_bitset<>, std::pair<double, boost::dynamic_bitset<>>> &coversets)
{
    std::unordered_map<boost::dynamic_bitset<>, std::vector<std::pair<boost::dynamic_bitset<>, double>>> containers_already_hit_by_unique_guide;

    log_buffer << "Monophonic threshold: " << monophonic_threshold << std::endl;

    if (monophonic_threshold > 0)
    {
        auto it = coversets.begin();
        while (it != coversets.end())
        {
            boost::dynamic_bitset<> new_guide = it->first;
            double new_score = it->second.first;
            boost::dynamic_bitset<> container_bitset = it->second.second;

            // Below, we want to keep only m guides per container where
            // that guide hits only this set of container and none other.
            // If a container (or set of containers) is already hit by an mp-hitting-guide and we encounter
            // another one, mark the new guide for deletion from coversets and skip adding it to hit_species.
            if (container_bitset.count() <= monophonic_threshold)
            {
                // and if this container already has m representative guides that hits it...
                auto container_bitset_iterator = containers_already_hit_by_unique_guide.find(container_bitset);
                if (container_bitset_iterator != containers_already_hit_by_unique_guide.end())
                {
                    if (container_bitset_iterator->second.size() == cut_multiplicity)
                    {
                        // Get the vector of representative guides and their scores in pairs
                        std::vector<std::pair<boost::dynamic_bitset<>, double>> represent_guides_vector = container_bitset_iterator->second;

                        std::pair<boost::dynamic_bitset<>, double> smallest_pair;
                        boost::dynamic_bitset<> smallest_old_guide = represent_guides_vector[0].first;
                        double smallest_score = represent_guides_vector[0].second;

                        // Find the weakest representative guide
                        for (const std::pair<boost::dynamic_bitset<>, double> &pair : represent_guides_vector)
                        {
                            if (pair.second < smallest_score)
                            {
                                smallest_pair = pair;
                                smallest_score = pair.second;
                                smallest_old_guide = pair.first;
                            }
                        }

                        // If the representative guide we saw earlier is a worse cutter, flag its score.
                        if (new_score > smallest_score)
                        {
                            auto it = coversets.find(smallest_old_guide);
                            it->second.first = 0; // A guide with a score of 0 will be removed in kirschtorte.cpp.
                                                // Removing a previous object here will mess with the for-loop that we are in.
                        }
                        else // Else if the new guide is not better than the weakest rep guide, discard the new guide and continue.
                        {
                            it = coversets.erase(it); // This erases the current element and returns an iterator to the next object.
                            continue; // Therefore do not update the iterator (it++) below.
                        }
                    }
                    else
                    {
                        // Vector is not full yet. Add the new guide to it.
                        // I.e., the multiplicity count is not met yet. This container still needs to be hit more times.
                        containers_already_hit_by_unique_guide[container_bitset].push_back(std::pair<boost::dynamic_bitset<>, double>(new_guide, new_score));
                    }
                }
                // If this species still needs representative guides...
                else
                {
                    // Indicate that this species has one representative guide now.
                    std::vector<std::pair<boost::dynamic_bitset<>, double>> vector_of_representatives;
                    vector_of_representatives.push_back(std::pair<boost::dynamic_bitset<>, double>(new_guide, new_score));
                    containers_already_hit_by_unique_guide[container_bitset] = vector_of_representatives;
                }
            }

            it++;
        }
    }
}

// Function to XOR every two bits and count mismatches
int xor_and_count_mismatches(const boost::dynamic_bitset<>& bitset1, const boost::dynamic_bitset<>& bitset2)
{
    int mismatches = 0;

    for (size_t i = 0; i < bitset1.size(); i += 2) {
        // XOR each pair of bits and check if the result is 1 (mismatch)
        bool result = (bitset1[i] != bitset2[i]) || (bitset1[i + 1] != bitset2[i + 1]);
        
        if (result) {
            ++mismatches;
        }
    }
    return mismatches;
}


// void decorate_with_clustering_median(
//     std::size_t seed_length,
//     std::size_t cut_multiplicity,  // aka m guides required per container
//     std::size_t mismatched_allowed_after_seed,
//     std::ostringstream &log_buffer,
//     std::map<boost::dynamic_bitset<>, std::pair<double, boost::dynamic_bitset<>>> &coversets)
// {
//     std::size_t left_side_bits_length;
//     std::size_t seed_bits_length = seed_length * 2;  // Each nucleotide is represented with 2 bits.
//     std::map<boost::dynamic_bitset<>, std::set<boost::dynamic_bitset<>>> seed_to_set;

//     auto it = coversets.begin();
//     while (it != coversets.end())
//     {
//         double new_score = it->second.first;

//         // Guide already marked for removal. Don't bother.
//         if (new_score <= 0)
//         {
//             it++;
//             continue;
//         }

//         boost::dynamic_bitset<> new_guide = it->first;

//         //  {'A', "00"}, {'G', "01"}, {'C', "10"}, {'T', "11"}
//         // AAAAAAGCTTGCTCTTTGCC 
//         // 0000000000000110111101101110111111011010
//         //
//         // Left side bits -        Seed bits
//         // AAAAAAGC         T T G C T C T T T G C C
//         // 0000000000000110 111101101110111111011110
//         left_side_bits_length = new_guide.size() - seed_bits_length;
//         boost::dynamic_bitset<> seed_bits(seed_bits_length, new_guide.to_ulong() & ((1 << seed_bits_length) - 1));
//         boost::dynamic_bitset<> new_guide_left_side_bits(left_side_bits_length, new_guide.to_ulong() >> (new_guide.size() - left_side_bits_length));

//         const auto it_seed_to_set = seed_to_set.find(seed_bits);
//         if (it_seed_to_set != seed_to_set.end())
//         {
//             seed_to_set[seed_bits].insert(new_guide_left_side_bits);
//         }
//         // First in line.
//         else
//         {
//             seed_to_set[seed_bits] = std::set<boost::dynamic_bitset<>>{new_guide_left_side_bits};
//         }

//         it++;
//     }

//     for (const auto map_it : seed_to_set)
//     {
//         std::set<boost::dynamic_bitset<>> left_sides_with_same_seed = map_it->second;

//         // Find the median string
//         // Variable to store bitsets with mismatched_allowed_after_seed or fewer mismatches to each other bitset
//         std::vector<boost::dynamic_bitset<>> candidates;
//         boost::dynamic_bitset<> left_side_of_candidate_guide(new_guide_left_side_bits.size());

//         bool found_a_median = false;

//         // Iterate through each bitset and compare it with every other bitset
//         for (const auto& bitset_outer : left_sides_with_same_seed)
//         {
//             bool valid_candidate = true;

//             for (const auto& bitset_inner : left_sides_with_same_seed)
//             {
//                 if ((bitset_outer != bitset_inner) && (xor_and_count_mismatches(bitset_outer, bitset_inner) > mismatched_allowed_after_seed))
//                 {
//                     valid_candidate = false;
//                     break;
//                 }
//             }
//             if (valid_candidate)
//             {
//                 found_a_median = true;
//                 left_side_of_candidate_guide = bitset_outer;
//             }
//         }

//         // if found a median
//         // else if not found a median
        
//     }

// }

void decorate_with_clustering(
    std::size_t seed_length,
    std::size_t cut_multiplicity,  // aka m guides required per container
    std::size_t mismatched_allowed_after_seed,
    std::ostringstream &log_buffer,
    std::map<boost::dynamic_bitset<>, std::pair<double, boost::dynamic_bitset<>>> &coversets)
{
    std::size_t left_side_bits_length;
    std::size_t seed_bits_length = seed_length * 2;  // Each nucleotide is represented with 2 bits.
    std::map<boost::dynamic_bitset<>, std::vector<boost::dynamic_bitset<>>> seed_to_vec;

    auto it = coversets.begin();
    while (it != coversets.end())
    {
        double new_score = it->second.first;

        // Guide already marked for removal. Don't bother.
        if (new_score <= 0)
        {
            it++;
            continue;
        }

        boost::dynamic_bitset<> new_guide = it->first;
        boost::dynamic_bitset<> container_bitset = it->second.second;

        //  {'A', "00"}, {'G', "01"}, {'C', "10"}, {'T', "11"}
        // AAAAAAGCTTGCTCTTTGCC 
        // 0000000000000110111101101110111111011010
        //
        // Left side bits -        Seed bits
        // AAAAAAGC        T T G C T C T T T G C C
        // 0000000000000110111101101110111111011110
        left_side_bits_length = new_guide.size() - seed_bits_length;
        boost::dynamic_bitset<> seed_bits(seed_bits_length, new_guide.to_ulong() & ((1 << seed_bits_length) - 1));
        boost::dynamic_bitset<> new_guide_left_side_bits(left_side_bits_length, new_guide.to_ulong() >> (new_guide.size() - left_side_bits_length));

        auto it_seed_to_vec = seed_to_vec.find(seed_bits);
        if (it_seed_to_vec != seed_to_vec.end())
        {
            bool found_an_heir = false;
            std::size_t fewest_mismatches = mismatched_allowed_after_seed + 1;
            boost::dynamic_bitset<> reigning_guide_w_fewest_mismatches(new_guide.size());
            boost::dynamic_bitset<> left_side_of_reigning_guide_w_fewest_mismatches(left_side_bits_length);
            std::vector<boost::dynamic_bitset<>> left_sides_with_same_seed = it_seed_to_vec->second;

            // Iterate over the vector in seed_to_vec and find the sequence
            // with the fewest mismatches within accepted threshold (mismatched_allowed_after_seed).
            for (boost::dynamic_bitset<>& reigning_guide_left_side : left_sides_with_same_seed)
            {
                // Perform XOR on every two bits and count the mismatches.
                std::size_t mismatches_after_seed = xor_and_count_mismatches(new_guide_left_side_bits, reigning_guide_left_side);

                if ((mismatches_after_seed <= mismatched_allowed_after_seed) && (mismatches_after_seed < fewest_mismatches))
                {
                    found_an_heir = true;

                    // New fewest matches
                    fewest_mismatches = mismatches_after_seed;
                    left_side_of_reigning_guide_w_fewest_mismatches = reigning_guide_left_side; 
                }
            }

            // None of the reigning guides have few enough mismatches to the new guide.
            if (found_an_heir == false)
            {
                // The new guide ('s left side bits) is now also a reigning.
                seed_to_vec[seed_bits].push_back(new_guide_left_side_bits);
            }
            // Get the best reigning guide and its targets.
            else
            {
                // Mark the new guide for deletion.
                coversets[new_guide].first = 0;

                //
                // Reconstruct the reigning guide
                //
                
                // Copy left side of reigning.
                boost::dynamic_bitset<> reigning_guide_left_side_extended(left_side_of_reigning_guide_w_fewest_mismatches);
                // Extend it.
                reigning_guide_left_side_extended.resize(new_guide.size());
                // Make room for seed_bits.
                reigning_guide_left_side_extended <<= seed_bits.size();

                // Copy seed_bits.
                boost::dynamic_bitset<> seed_bits_extended(seed_bits);
                // Extend it.
                seed_bits_extended.resize(reigning_guide_left_side_extended.size());

                // Reconstruct the reigning guide by concatenating left and seed.
                reigning_guide_w_fewest_mismatches = reigning_guide_left_side_extended | seed_bits_extended;

                // The reigning guide inherits the targets of the new guide.
                coversets[reigning_guide_w_fewest_mismatches].second |= coversets[new_guide].second;
            }
        }
        // First reigning.
        else
        {
            seed_to_vec[seed_bits] = std::vector<boost::dynamic_bitset<>>{new_guide_left_side_bits};
        }

        it++;
    }
}

void decorate_with_clustering_multiplicity(
    std::size_t seed_length,
    std::size_t cut_multiplicity,  // aka m guides required per container
    std::size_t mismatched_allowed_after_seed,
    std::ostringstream &log_buffer,
    std::map<boost::dynamic_bitset<>, std::pair<double, boost::dynamic_bitset<>>> &coversets)
{

    std::size_t left_side_bits_length;
    std::size_t seed_bits_length = seed_length * 2;  // Each nucleotide is represented by 2 bits.

    std::map<boost::dynamic_bitset<>, std::vector<boost::dynamic_bitset<>>> seed_to_vec;

    auto it = coversets.begin();
    while (it != coversets.end())
    {
        double new_score = it->second.first;

        if (new_score <= 0)
        {
            it++;
            continue;
        }

        boost::dynamic_bitset<> new_guide = it->first;
        boost::dynamic_bitset<> container_bitset = it->second.second;

        //  {'A', "00"}, {'G', "01"}, {'C', "10"}, {'T', "11"}
        // AAAAAAGCTTGCTCTTTGCC 
        // 0000000000000110111101101110111111011010
        //
        // Left side bits -        Seed bits
        // AAAAAAGC        T T G C T C T T T G C C
        // 0000000000000110111101101110111111011110
        std::size_t left_side_bits_length = new_guide.size() - seed_bits_length;
        boost::dynamic_bitset<> seed_bits(seed_bits_length, new_guide.to_ulong() & ((1 << seed_bits_length) - 1));
        boost::dynamic_bitset<> new_guide_left_side_bits(left_side_bits_length, new_guide.to_ulong() >> (new_guide.size() - left_side_bits_length));

        auto it_seed_to_vec = seed_to_vec.find(seed_bits);
        if (it_seed_to_vec != seed_to_vec.end())
        {
            bool found_an_heir = false;
            std::vector<boost::dynamic_bitset<>> left_sides_with_same_seed = it_seed_to_vec->second;

            // Iterate over the vector in seed_to_vec and find the sequence.
            // with the fewest mismatches within accepted threshold (mismatched_allowed_after_seed).
            for (boost::dynamic_bitset<>& reigning_guide_left_side : left_sides_with_same_seed)
            {
                // Perform XOR on every two bits and count the mismatches.
                std::size_t mismatches_after_seed = xor_and_count_mismatches(new_guide_left_side_bits, reigning_guide_left_side);

                if (mismatches_after_seed <= mismatched_allowed_after_seed)
                {
                    found_an_heir = true;

                    // Copy reigning
                    boost::dynamic_bitset<> reigning_guide_left_side_extended(reigning_guide_left_side);
                    // Extend it
                    reigning_guide_left_side_extended.resize(reigning_guide_left_side.size() + seed_bits.size());
                    // Make room for seed_bits
                    reigning_guide_left_side_extended <<= seed_bits.size();

                    // Copy seed_bits
                    boost::dynamic_bitset<> seed_bits_extended(seed_bits);
                    // Extend it
                    seed_bits_extended.resize(reigning_guide_left_side_extended.size());

                    // Reconstruct the reigning guide
                    boost::dynamic_bitset<> reigning_guide_w_fewest_mismatches(reigning_guide_left_side_extended | seed_bits_extended); 

                    // The reigning guide inherits the targets of the new guide.
                    coversets[reigning_guide_w_fewest_mismatches].second |= coversets[new_guide].second;
                }
            }

            // None of the reigning guides have few enough mismatches to the new guide.
            if (found_an_heir == false)
            {
                // The new guide is also an reigning.
                seed_to_vec[seed_bits].push_back(new_guide_left_side_bits);
            }
            else
            {
                // Mark the new guide for deletion.
                coversets[new_guide].first = 0;
            }
        }
        // First reigning.
        else
        {
            // std::vector<boost::dynamic_bitset<>> vec;
            // vec.push_back(new_guide_left_side_bits);
            seed_to_vec[seed_bits] = std::vector<boost::dynamic_bitset<>>{new_guide_left_side_bits};
        }

        it++;
    }
}