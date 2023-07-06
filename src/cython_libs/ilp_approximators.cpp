#include "allegro/definitions.h"
#include "allegro/decode_bitset.h"
#include "allegro/ilp_approximators.h"

std::vector<GuideStruct> randomized_rounding(
    std::vector<operations_research::MPVariable *> &feasible_solutions,
    boost::dynamic_bitset<> all_containers_bitset,
    std::map<boost::dynamic_bitset<>, std::pair<double, boost::dynamic_bitset<>>> &coversets,
    std::size_t multiplicity,
    std::size_t num_trials,
    std::ostringstream &log_buffer)
{
    std::cout << BLUE << "> " << RESET << "Using randomized rounding with " << num_trials << " trials." << std::endl;
    log_buffer << "Using randomized rounding with " << num_trials << " trials." << std::endl;

    std::size_t original_size = all_containers_bitset.size();
    std::size_t new_size = original_size * multiplicity;

    boost::dynamic_bitset<> extended_all_containers_bitset(new_size);
    all_containers_bitset.resize(new_size);

    for (std::size_t i = 0; i < multiplicity; i++)
    {
        extended_all_containers_bitset <<= original_size;
        extended_all_containers_bitset |= all_containers_bitset;
    }

    std::set<std::string> winners;
    std::size_t len_winners = INT_MAX; // Inifinity, initially
    std::size_t trial_with_smallest_size = 0;

    // On each invocation, dist(rng) returns a random floating-point value
    // uniformly distributed in the range [min, max).
    boost::random::mt19937 rng(std::time(nullptr));
    boost::random::uniform_real_distribution<double> dist(0, 1);

    std::size_t trial = 1;
    while (trial < num_trials + 1)
    {
        if (trial % 1000 == 0)
        {
            // Writes on the same line
            std::cout << BLUE << "\r> " << RESET << "Trial " << trial << std::flush;
        }

        // Guide containers to cover
        boost::dynamic_bitset<> I_this_trial(new_size);
        std::set<std::string> winners_this_trial;
        std::size_t iterations_this_trial = 100000; // while-loop exit condition in case of bad luck

        while ((I_this_trial != extended_all_containers_bitset) && (iterations_this_trial != 0))
        {
            iterations_this_trial--;

            for (auto var_ptr : feasible_solutions)
            {
                operations_research::MPVariable *var = var_ptr;

                if (winners_this_trial.find(var->name()) != winners_this_trial.end())
                {
                    continue;
                }

                if ((var->solution_value() == 1.0) || (var->solution_value() > dist(rng)))
                {
                    // Encoded binary DNA sequence
                    boost::dynamic_bitset<> bitset(var->name());

                    // The species bit vector hit by this binary DNA sequence
                    boost::dynamic_bitset<> species_hit_by_this_guide = coversets[bitset].second;

                    species_hit_by_this_guide.resize(new_size);
                    boost::dynamic_bitset<> old_I(new_size);

                    while (species_hit_by_this_guide.find_first() != boost::dynamic_bitset<>::npos)
                    {
                        old_I = I_this_trial;
                        I_this_trial |= species_hit_by_this_guide;
                        species_hit_by_this_guide &= old_I;
                        species_hit_by_this_guide <<= original_size;
                    }

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

        trial += 1;
    }

    if (trial > 1000)
    {
        std::cout << std::endl;
    }

    std::vector<GuideStruct> decoded_winners;

    std::cout <<  BLUE << "> " << RESET << "Winners are:" << std::endl;
    log_buffer << "Winners are:" << std::endl;
    for (auto winner_str : winners)
    {
        boost::dynamic_bitset<> bitset(winner_str);

        double score = coversets[bitset].first;
        boost::dynamic_bitset<> species_hit_by_this_guide = coversets[bitset].second;

        std::string buffer;
        boost::to_string(species_hit_by_this_guide, buffer);

        std::string decoded_bitset = decode_bitset(winner_str);

        GuideStruct guide;
        guide.sequence = decoded_bitset;
        guide.score = score;
        guide.species_hit = buffer;

        decoded_winners.push_back(guide);

        std::cout << BLUE << "> " << RESET << decoded_bitset << std::endl;
        log_buffer << decoded_bitset << std::endl;
    }

    return decoded_winners;
}