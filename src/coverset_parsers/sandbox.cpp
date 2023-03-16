#include <iostream>
#include <string>
#include <unordered_set>
#include <boost/random.hpp>

int main() {
  // Create a random number generator and distribution
  boost::random::mt19937 rng(std::time(nullptr));
  boost::random::uniform_real_distribution<> dist(0, 1);

  // Generate a random number between 0 and 1
  float random = dist(rng);

  // Print the random number to the console
  std::cout << "The random number is: " << random << std::endl;

  std::unordered_set<std::string> empty;
  empty.insert("NA");

  std::cout << *empty.begin() << std::endl;

  return 0;
}