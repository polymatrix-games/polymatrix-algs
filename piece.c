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

    matrix_t *xt, *xht, *vi_x, *vi_xh, *tmp;
    k = 0;
    //strategy_print(x);
    //strategy_print(xp);
    //strategy_print(xh);
    for (i = 0; i < x->players; ++i){
        xt = matrix_trans(x->strategies[i]);
        xht = matrix_trans(xh->strategies[i]);
        //printf("xh\n");
        //matrix_print(xht);
        vi_x = polymatrix_compute_vi(game, x, i);
        vi_xh = polymatrix_compute_vi(game, xh, i);

        tmp = matrix_mul_vec_mat(xt, vi_x);
        c2 = tmp->data[0][0];
        matrix_free(tmp);

        tmp = matrix_mul_vec_mat(xt, vi_xh);
        //printf("vi_xh\n");
        //matrix_print(vi_xh);
        e2 = tmp->data[0][0];
        matrix_free(tmp);

        tmp = matrix_mul_vec_mat(xht, vi_x);
        //printf("vi_x\n");
        //matrix_print(vi_x);
        e3 = tmp->data[0][0];
        matrix_free(tmp);

        tmp = matrix_mul_vec_mat(xht, vi_xh);
        e_2 = tmp->data[0][0];
        matrix_free(tmp);

        //printf("Piece\n============\n");
        for (l = 0; l < scount[i]; ++l, ++k) {
            p->data[k][0] = vi_x->data[l][0] - c2;
            //printf("(%lf - %lf) ", vi_x->data[l][0], c2);
            p->data[k][1] = vi_xh->data[l][0] - e2 - e3;
            //printf("+ eps (%lf - %lf - %lf) ", vi_xh->data[l][0], e2 , e3);
            p->data[k][2] = -e_2;
            //printf("+ eps2 %lf \n", e_2);
        }
        //printf("\n============\n");

        matrix_free(vi_x);
        matrix_free(vi_xh);
    }
    //matrix_print(p);

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
    fprintf(f, "plot max(");
    for (i = 0; i < piece->nrows; ++i) {
        for (j = 0; j < piece->ncols; ++j) {
            fprintf(f, " + %lf*x**%d", piece->data[i][j], j);
        }
        fprintf(f,", ");
    }
    fprintf(f, "0)\npause -1");
}
