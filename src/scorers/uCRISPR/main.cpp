#include "include/uCRISPR_scorer.h"

#include <string>
#include <vector>
#include <iostream>

int main(int argc, char* argv[]) {
    std::string input;
    std::vector<std::string> sequences_vec;

    while (std::getline(std::cin, input)) {
        sequences_vec.push_back(input);
    }

    uCRISPR_scorer::uCRISPR_scorer scorer;
    std::vector<double> scores_vec = scorer.score_guides(sequences_vec);

    for (int i = 0; i < sequences_vec.size(); i++) {
        std::cout << sequences_vec[i] + " " + std::to_string(scores_vec[i]) << std::endl;
    }

    return 0;
}
