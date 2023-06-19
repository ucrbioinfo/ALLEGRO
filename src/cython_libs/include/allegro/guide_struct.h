#ifndef GUIDESTRUCT_H
#define GUIDESTRUCT_H

#include <string>

struct GuideStruct {
    std::string sequence;
    double score;
    std::string species_hit; // looks like this "010101001"
};

#endif