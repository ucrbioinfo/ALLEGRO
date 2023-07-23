#include "include/uCRISPR_scorer.h"

#include <map>
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
    scorer.score_guides(sequences_vec);

    return 0;
}
