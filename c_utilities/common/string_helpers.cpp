#include "string_helpers.h"
#include <string.h>

std::vector<int> split_string(std::string input, const char *delimiter)
{
	std::vector<int> ret_values;
	char *pch;
	char *str = strdup(input.c_str());
	pch = strtok(str, delimiter);
	while (pch != NULL) {
		int ind = strtod(pch, NULL);
		ret_values.push_back(ind);
		pch = strtok(NULL, delimiter);
	}
	free(str);
	return ret_values;
}

std::string remove_file_ending(std::string filename); //from readimage.cpp

std::string get_basename(std::string input)
{
	input = remove_file_ending(input);
	char *input_copy = strdup(input.c_str());
	char *basename_str = basename(input_copy);
	std::string output_basename = std::string(basename_str);
	free(input_copy);
	return output_basename;
}
