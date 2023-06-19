#ifndef UCRISPR_SCORER_H
#define UCRISPR_SCORER_H

#include <map>
#include <vector>

namespace uCRISPR_scorer
{
    class uCRISPR_scorer
    {
    public:
        uCRISPR_scorer();
        ~uCRISPR_scorer();

        // Call this function with a vector of length 23 guide RNA sequences including the cas9 PAM
        // E.g., AGCGTACCCCCAGGTCTTGCAGG
        std::vector<double> score_guides(std::vector<std::string> sequences);
        double score_guide(std::string sequence);

    private:
        std::map<std::string, double> parameters_on;
        
        void ReadParameters();
        double GetSgRNAEnergy(std::string sequence);
    };
}

#endif