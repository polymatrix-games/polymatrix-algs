#include "../util/matrix.h"
#include "../util/polymatrix.h"
#include "../util/strategy.h"

#ifndef PIECE_H
#define PIECE_H

matrix_t *piece_create(polymatrix_t *game, strategy_profile_t *x, strategy_profile_t *xp, int *scount);
matrix_t *piece_p(matrix_t *piece);
void print_piece(matrix_t *piece, FILE *f);
double piece_eval(matrix_t *p, double eps);

#endif
