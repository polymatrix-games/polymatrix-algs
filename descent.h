#include "../util/polymatrix.h"

#ifndef DESCENT_H
#define DESCENT_H

extern int line_search;
void solve_descent(polymatrix_t *game, double delta);
strategy_profile_t *run_descent(polymatrix_t *game, double delta);

#endif
