#ifndef PROGRESS_BAR_H_DEFINED
#define PROGRESS_BAR_H_DEFINED

#include <string>

/**
 * Print progress bar.
 *
 * \param action String describing the action taken
 * \param last_progress Last progress sent to the function
 * \param progress Current progress
 * \param total Total progress that can be reached
 **/
void print_progress_bar(std::string action, int last_progress, int progress, int total);

#endif
