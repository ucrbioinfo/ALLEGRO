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

        double score_guide(const std::string &sequence);
        void score_guides(const std::vector<std::string> &sequences);

    private:
        std::map<std::string, double> parameters_on;

        void ReadParameters();
        
        double GetSgRNAEnergy(const std::string &sequence);
        void threaded_scorer(const std::vector<std::string>& sequences, std::vector<double>& scores, size_t start, size_t end);
    };
}

#endif