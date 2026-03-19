#include "allegro/definitions.h"
#include "allegro/decode_bitset.h"

std::string decode_bitset(boost::dynamic_bitset<> encoded)
{
    std::string decoded = "";
    std::string buffer;
    std::string next_two_bases;

    boost::to_string(encoded, buffer);

    for (std::size_t i = 0; i < buffer.length(); i += 2)
    {
        next_two_bases = buffer.substr(i, 2);

        if (next_two_bases == sA_SHIFT)
        {
            decoded += "A";
        }
        else if (next_two_bases == sC_SHIFT)
        {
            decoded += "C";
        }
        else if (next_two_bases == sG_SHIFT)
        {
            decoded += "G";
        }
        else if (next_two_bases == sT_SHIFT)
        {
            decoded += "T";
        }
    }

    return decoded;
}

std::string decode_bitset(const std::string &encoded_str)
{
    std::string decoded = "";
    std::string next_two_bases;

    for (std::size_t i = 0; i < encoded_str.length(); i += 2)
    {
        next_two_bases = encoded_str.substr(i, 2);

        if (next_two_bases == sA_SHIFT)
        {
            decoded += "A";
        }
        else if (next_two_bases == sC_SHIFT)
        {
            decoded += "C";
        }
        else if (next_two_bases == sG_SHIFT)
        {
            decoded += "G";
        }
        else if (next_two_bases == sT_SHIFT)
        {
            decoded += "T";
        }
    }

    return decoded;
}