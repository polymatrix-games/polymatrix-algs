#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "piece.h"

matrix_t *piece_create(polymatrix_t *game, strategy_profile_t *x, strategy_profile_t *xp, int *scount)
{
    int i, j, k, l, total = 0;
    for (i = 0; i < game->players; ++i)
        total += scount[i];

    matrix_t *p = matrix_alloc(total, 3);

    strategy_profile_t *xh = strategy_sub(xp, x, scount);

    double c2, e2, e3, e_2;

    matrix_t *xt, *xht,*Aij_xj, *Aij_xhj, *tmp;
    k = 0;
    for (i = 0; i < x->players; ++i){
        xt = matrix_trans(x->strategies[i]);
        xht = matrix_trans(xh->strategies[i]);
        Aij_xj = matrix_alloc(scount[i], 1);
        Aij_xhj = matrix_alloc(scount[i], 1);
        for (j = 0; j < x->players; ++j) {
            if (!game->payoffs[i][j])
                continue;
            tmp = matrix_mul(game->payoffs[i][j], x->strategies[j]);
            matrix_add_dst(Aij_xj, tmp, Aij_xj);
            matrix_free(tmp);

            tmp = matrix_mul(game->payoffs[i][j], xh->strategies[j]);
            matrix_add_dst(Aij_xhj, tmp, Aij_xhj);
            matrix_free(tmp);
        }
        tmp = matrix_mul_vec_mat(xt, Aij_xj);
        c2 = tmp->data[0][0];
        matrix_free(tmp);

        tmp = matrix_mul_vec_mat(xt, Aij_xhj);
        e2 = tmp->data[0][0];
        matrix_free(tmp);

        tmp = matrix_mul_vec_mat(xht, Aij_xj);
        e3 = tmp->data[0][0];
        matrix_free(tmp);

        tmp = matrix_mul_vec_mat(xht, Aij_xhj);
        e_2 = tmp->data[0][0];
        matrix_free(tmp);

        for (l = 0; l < scount[i]; ++l, ++k) {
            p->data[k][0] = Aij_xj->data[l][0] - c2;
            p->data[k][1] = Aij_xhj->data[l][0] - e2 - e3;
            p->data[k][2] = e_2;
        }

        matrix_free(Aij_xj);
        matrix_free(Aij_xhj);
    }

    return p;
}

matrix_t *piece_p(matrix_t *piece)
{
    matrix_t *p = matrix_alloc(piece->nrows, piece->ncols - 1);

    int i, j;
    for (i = 0; i < p->nrows; ++i) {
        for (j = 0; j < p->ncols; ++j)
            p->data[i][j] = (j + 1) * piece->data[i][j+1];
    }

    return p;
}

double piece_eval(matrix_t *p, double eps)
{
    double m = 0;
    int i, j;

    for (i = 0; i < p->nrows; ++i){
        double t = 0;
        for (j = 0; j < p->ncols; ++j)
            t += p->data[i][j] * pow(eps, j);
        m = (t > m) ? t : m;
    }

    return m;
}

void print_piece(matrix_t *piece, FILE *f)
{
    int i, j;
    fprintf(f, "set nokey\nset xrange [0:1]\n");
    fprintf(f, "plot");
    for (i = 0; i < piece->nrows; ++i) {
        for (j = 0; j < piece->ncols; ++j) {
            fprintf(f, " + %lf*x**%d", piece->data[i][j], j);
        }
        fprintf(f,", ");
    }
    fprintf(f, "\npause -1");
}
