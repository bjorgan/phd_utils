#include "progress_bar.h"

void print_progress_bar(std::string action, int last_progress, int progress, int total) {
	if (progress == 0) {
		fprintf(stderr, "%s\n\n", action.c_str());
	}

	int MAX_NUM_SIGNS = 20;
	int num_signs = progress*1.0/(total*1.0)*MAX_NUM_SIGNS;
	int prev_num_signs = last_progress*1.0/(total*1.0)*MAX_NUM_SIGNS;

	if (prev_num_signs != num_signs) {
		char str[MAX_NUM_SIGNS+3] = {0};
		str[0] = '[';
		str[MAX_NUM_SIGNS] = ']';
		for (int i=1; i < MAX_NUM_SIGNS; i++) {
			if (i < num_signs) {
				str[i] = '#';
			} else {
				str[i] = ' ';
			}
		}
		fprintf(stderr, "\033[1A%s\n", str);
	}
}

