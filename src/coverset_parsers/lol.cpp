#include <iostream>
#include <boost/dynamic_bitset.hpp>


int main() {

    std::string seq = "ACTG";
    std::size_t seq_len = seq.length();

    boost::dynamic_bitset<> encoded(seq_len * 2);

    boost::dynamic_bitset<> C_SHIFT(seq_len * 2, 0b01);

    std::cout << "Bits: " << encoded << std::endl;


    encoded |= C_SHIFT;

    std::cout << "Bits: " << encoded << std::endl;

    std::string decoded = "";

    boost::dynamic_bitset<> LSBs(seq_len * 2, 0b11);
    LSBs &= encoded;

    std::string buffer;
    boost::to_string(LSBs, buffer);

    std::cout << buffer << std::endl;


    return 0;
}