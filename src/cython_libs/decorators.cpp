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

void decorate_with_clustering(
    std::ostringstream &log_buffer,
    std::map<boost::dynamic_bitset<>, std::pair<double, boost::dynamic_bitset<>>> &coversets)
{
    // TODO - preclustering
}