#ifndef STRING_HELPERS_H_DEFINED
#define STRING_HELPERS_H_DEFINED

#include <vector>
#include <string>

/**
 * Split input string at delimiter.
 *
 * \param input Input string
 * \param delimiter Delimiters
 * \return Split string
 **/
std::vector<int> split_string(std::string input, const char *delimiter);

/**
 * Remove file extension and paths.
 **/
std::string get_basename(std::string input);

#endif
