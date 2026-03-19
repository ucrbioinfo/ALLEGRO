#ifndef LOGGER_H
#define LOGGER_H

#include <sstream>
#include <fstream>
#include <iostream>
#include <filesystem>

void log_info(std::ostringstream &log_buffer, std::string output_directory);

#endif