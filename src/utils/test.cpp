#include <iostream>
#include <vector>
#include <string>
#include <thread>
#include <mutex>
#include <random>

std::mutex mtx; // Mutex to protect shared data

class lol {
public:
    // Function to calculate the score for a given sequence
    double score_sequence(const std::string& sequence) {
        // Simulating some score calculation
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<double> dis(0.0, 1.0);
        double score = dis(gen);

        return score;
    }

    // Callable object to invoke calculate_scores member function
    struct ScoreCalculator {
        lol* instance;
        const std::vector<std::string>& sequences;
        std::vector<double>& scores;
        size_t start;
        size_t end;

        void operator()() const {
            instance->calculate_scores(sequences, scores, start, end);
        }
    };

    // Thread function to calculate scores for a range of sequences
    void calculate_scores(const std::vector<std::string>& sequences, std::vector<double>& scores, size_t start, size_t end) {
        for (size_t i = start; i < end; ++i) {
            double score = score_sequence(sequences[i]);

            std::lock_guard<std::mutex> lock(mtx); // Lock the mutex to protect shared data
            scores[i] = score; // Store the score in the corresponding index
        }
    }

    void score() {
        std::vector<std::string> sequences = { "AAA", "BBB", "CCC", "DDD", "EEE" };
        std::vector<double> scores(sequences.size());

        // Number of threads to use
        size_t num_threads = std::thread::hardware_concurrency();

        // Create vector of threads
        std::vector<std::thread> threads;
        threads.reserve(num_threads);

        // Calculate chunk size for each thread
        size_t chunk_size = sequences.size() / num_threads;

        // Distribute the work among threads
        for (size_t i = 0; i < num_threads; ++i) {
            size_t start = i * chunk_size;
            size_t end = (i == num_threads - 1) ? sequences.size() : start + chunk_size;

            ScoreCalculator calculator{ this, sequences, scores, start, end };
            threads.emplace_back(calculator);
        }

        // Wait for all threads to finish
        for (auto& thread : threads) {
            thread.join();
        }

        // Print the sequences along with their scores
        for (size_t i = 0; i < sequences.size(); ++i) {
            std::cout << "Sequence: " << sequences[i] << ", Score: " << scores[i] << std::endl;
        }
    }
};

int main() {
    lol l;
    l.score();
    return 0;
}