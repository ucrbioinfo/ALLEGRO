#include "allegro/logger.h"


void log_info(std::ostringstream &log_buffer, std::string output_directory)
{
    std::filesystem::path dirpath(output_directory);
    std::filesystem::path filepath = dirpath / "solver_log.txt";

    std::ofstream log_file(filepath);

    if (!log_file.is_open())
    {
        std::cerr << "Unable to open log file!" << std::endl;
    }

    log_file << log_buffer.str();
    log_file.close();
}