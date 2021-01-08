#ifndef GETOPT_LONG_HELPERS_H_DEFINED
#define GETOPT_LONG_HELPERS_H_DEFINED

#include <getopt.h>

/**
 * Print program usage instructions to stdout.
 *
 * \param General usage instructions
 * \param long_options List of long options used in getopts_long
 * \param short_options List of short options used in getopts_long
 * \param option_descriptions Description of each option
 **/
void getopt_long_show_help(const char *usage_instructions, struct option long_options[], const char *short_options, const char **option_descriptions);

#endif
