#include "allegro/ilp_approximators.h"

std::vector<std::pair<std::string, std::string>> randomized_rounding(
    std::vector<operations_research::MPVariable *> &feasible_solutions,
    boost::dynamic_bitset<> &all_species_bitset,
    std::unordered_map<boost::dynamic_bitset<>, std::pair<char, boost::dynamic_bitset<>>> &coversets,
    std::size_t num_containers,
    std::size_t num_trials,
    std::ostringstream &log_buffer)
{
    // Algorithm source:
    // https://web.archive.org/web/20230325165759/https://theory.stanford.edu/~trevisan/cs261/lecture08.pdf
    std::cout << "Using randomized rounding with " << num_trials << " trials." << std::endl;
    log_buffer << "Using randomized rounding with " << num_trials << " trials." << std::endl;

    std::unordered_set<std::string> winners;
    std::size_t len_winners = INT_MAX; // Inifinity, initially
    std::size_t trial_with_smallest_size = 0;

    // On each invocation, dist(rng) returns a random floating-point value
    //  uniformly distributed in the range [min, max).
    boost::random::mt19937 rng(std::time(nullptr));
    boost::random::uniform_real_distribution<double> dist(0, 1);

    for (std::size_t trial = 1; trial < num_trials + 1; trial++)
    {
        // Guide containers to cover
        boost::dynamic_bitset<> I_this_trial(num_containers);
        std::unordered_set<std::string> winners_this_trial;
        std::size_t iterations_this_trial = 100000; // while-loop exit condition in case of bad luck

        while ((I_this_trial != all_species_bitset) && (iterations_this_trial != 0))
        {
            iterations_this_trial--;

            for (auto var_ptr : feasible_solutions)
            {
                operations_research::MPVariable *var = var_ptr;

                if ((var->solution_value() == 1.0) || (var->solution_value() > dist(rng)))
                {
                    // Encoded binary DNA sequence
                    boost::dynamic_bitset<> bitset(var->name());

                    // The species bit vector hit by this binary DNA sequence
                    boost::dynamic_bitset<> species_hit_by_this_guide = coversets[bitset].second;

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
    log_buffer << "Winners are:" << std::endl;
    for (auto winner_str : winners)
    {
        boost::dynamic_bitset<> bitset(winner_str);
        boost::dynamic_bitset<> species_hit_by_this_guide = coversets[bitset].second;

        std::string buffer;
        boost::to_string(species_hit_by_this_guide, buffer);

        std::string decoded_bitset = decode_bitset(winner_str);
        decoded_winners.push_back(std::pair<std::string, std::string>(decoded_bitset, buffer));

        std::cout << decoded_bitset << std::endl;
        log_buffer << decoded_bitset << std::endl;
    }

    return decoded_winners;
}