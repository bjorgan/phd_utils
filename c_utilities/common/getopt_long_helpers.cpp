#include "getopt_long_helpers.h"
#include <string>
#include <cstring>

/**
 * Returns true if specified option's value is a char within the short options char array (and has thus a short-hand version of the long option)
 *
 * \param short_options Array over short option chars
 * \param long_option Option struct to check
 * \returns True if option has a short option, false otherwise
 **/
bool has_short_option(const char *short_options, struct option long_option)
{
	const char *ptr = strchr(short_options, long_option.val);
	if (ptr == NULL) {
		return false;
	}
	return true;
}

const int MAX_NUM_CHARS_IN_OPTION_PRINTING = 255;

void getopt_long_show_help(const char *usage_instructions, struct option long_options[], const char *short_options, const char **option_descriptions)
{
	//find maximum length of long options
	int max_length = 0;
	int index = 0;
	while (true) {
		const char *name = long_options[index].name;
		if (name == 0) {
			break;
		}

		int length = strlen(name);
		if (length > max_length) {
			max_length = length;
		}
		index++;
	}
	max_length += 6; //extra signs that are displayed after the long option in --help
	max_length += 8;

	//display initial description
	index = 0;
	int description_index = 0;
	printf("%s\n", usage_instructions);
	while (true) {
		char option_print_full[MAX_NUM_CHARS_IN_OPTION_PRINTING] = {0};
		std::fill(option_print_full, option_print_full + MAX_NUM_CHARS_IN_OPTION_PRINTING, ' ');
		char *option_print = option_print_full;

		if (long_options[index].name == 0) {
			break;
		}

		//display short option
		if (has_short_option(short_options, long_options[index])) {
			option_print += snprintf(option_print, MAX_NUM_CHARS_IN_OPTION_PRINTING, " -%c,", long_options[index].val);
		} else {
			option_print += snprintf(option_print, MAX_NUM_CHARS_IN_OPTION_PRINTING, "    ", long_options[index].val);
		}

		//display long option
		option_print += snprintf(option_print, MAX_NUM_CHARS_IN_OPTION_PRINTING, "--%s", long_options[index].name);

		//display argument
		switch (long_options[index].has_arg) {
			case required_argument:
				option_print += snprintf(option_print, MAX_NUM_CHARS_IN_OPTION_PRINTING, "=ARG");
			break;

			case optional_argument:
				option_print += snprintf(option_print, MAX_NUM_CHARS_IN_OPTION_PRINTING, "=OPT_ARG");
			break;
		}
		*option_print = ' ';

		option_print = option_print_full + max_length;
		if (option_descriptions[description_index] != NULL) {
			snprintf(option_print, MAX_NUM_CHARS_IN_OPTION_PRINTING-max_length, " %s.", option_descriptions[description_index]);
			description_index++;
		} else {
			snprintf(option_print, MAX_NUM_CHARS_IN_OPTION_PRINTING-max_length, " Documentation missing.");
		}
		printf("%s\n", option_print_full);
		index++;
	}
}
