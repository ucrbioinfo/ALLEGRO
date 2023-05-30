#include "allegro/decorators.h"

void decorate_with_monophonic(
    std::size_t cut_multiplicity,
    std::size_t monophonic_threshold,
    std::ostringstream &log_buffer,
    std::unordered_map<boost::dynamic_bitset<>, std::pair<char, boost::dynamic_bitset<>>> &coversets)
{
    std::unordered_map<boost::dynamic_bitset<>, std::vector<std::pair<boost::dynamic_bitset<>, char>>> containers_already_hit_by_unique_guide;

    log_buffer << "Monophonic threshold: " << monophonic_threshold << std::endl;

    if (monophonic_threshold > 0)
    {
        auto it = coversets.begin();
        while (it != coversets.end())
        {
            boost::dynamic_bitset<> new_guide = it->first;
            unsigned char new_score = it->second.first;
            boost::dynamic_bitset<> container_bitset = it->second.second;

            // Below, we want to keep only m guides per species where that guide hits
            // only this set of species and none other.
            // If a species (or set of species) is already hit by a one-hitting-guide and we encounter
            // another one, mark the new guide for deletion from coversets and
            // skip adding it to hit_species.
            if (container_bitset.count() <= monophonic_threshold)
            {
                // and if this species already has #multiplicity representative guides that hits it...
                auto container_bitset_iterator = containers_already_hit_by_unique_guide.find(container_bitset);
                if (container_bitset_iterator != containers_already_hit_by_unique_guide.end())
                {
                    if (container_bitset_iterator->second.size() == cut_multiplicity)
                    {
                        // Get the vector of representative guides and their scores in pairs
                        std::vector<std::pair<boost::dynamic_bitset<>, char>> represent_guides_vector = container_bitset_iterator->second;

                        std::pair<boost::dynamic_bitset<>, char> smallest_pair;
                        boost::dynamic_bitset<> smallest_old_guide = represent_guides_vector[0].first;
                        char smallest_score = represent_guides_vector[0].second;

                        for (const std::pair<boost::dynamic_bitset<>, char> &pair : represent_guides_vector)
                        {
                            if (pair.second < smallest_score)
                            {
                                smallest_pair = pair;
                                smallest_score = pair.second;
                                smallest_old_guide = pair.first;
                            }
                        }

                        // if the representative guide we saw earlier is a worse cutter, flag its score.
                        if (new_score > smallest_score)
                        {
                            auto it = coversets.find(smallest_old_guide);
                            it->second.first = 0;
                        }
                        // Else if the new guide is not better, delete it and continue.
                        else
                        {
                            it = coversets.erase(it);
                            continue; // Do not update the iterator below. The line above does that instead.
                        }
                    }
                    else
                    {
                        // Vector is not full yet. Add the new guide to it.
                        // I.e., the multiplicity count is not met yet. This container still needs to be hit more times.
                        containers_already_hit_by_unique_guide[container_bitset].push_back(std::pair<boost::dynamic_bitset<>, char>(new_guide, new_score));
                    }
                }
                // If this species still needs a representative guide...
                else
                {
                    // Indicate that this species has a representative guide now.
                    std::vector<std::pair<boost::dynamic_bitset<>, char>> vector_of_representatives;
                    vector_of_representatives.push_back(std::pair<boost::dynamic_bitset<>, char>(new_guide, new_score));
                    containers_already_hit_by_unique_guide[container_bitset] = vector_of_representatives;
                }
            }

            it++;
        }
    }
}