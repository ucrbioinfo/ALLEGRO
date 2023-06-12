#include <boost/dynamic_bitset.hpp>
#include <iostream>

int main() {
    boost::dynamic_bitset<> bitset(5, 0b11111);
    std::cout << "Original bitset: " << bitset << std::endl;

    return 0;
}
