
//#define _POSIX_C_SOURCE 199309L
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <pthread.h>
#include <time.h>

#include "../util/matrix.h"
#include "descent.h"
#include "piece.h"
#include "../util/strategy.h"

static int lp_cols = 0;
#define TIMEOUT 10*60
static int pivcount;
static double g_eps;
static strategy_profile_t *gx;
int line_search = -1;

double piece_time = 0;

int get_K_x(double eps, double *regret, int n, int *K_x)
{
    int i, count;
    count = 0;

    for (i = 0; i < n; ++i){
        //printf("%d %lf %lf\n", i, regret[i], eps);
        //if (err_eq(regret[i], eps))
        K_x[count++] = i;
    }

    return count;
}

void construct_C1_i(polymatrix_t *game, int *scount, int i, int sKx, int k, matrix_t *dbr_i, matrix_t *C, int r_off)
{
    int j;
    for (j = 0; j < dbr_i->nrows; ++j)
        C->data[r_off + j][k + 1] = -1;
    
    //printf("Sdri %d\n", dbr_i->nrows);
    int c_off = sKx + 1;
    for (j = 0; j < game->players; ++j) {
        matrix_t *A_ij;
        if (!game->payoffs[i][j]){
            A_ij = matrix_alloc(dbr_i->nrows, scount[j]);
        } else {
            A_ij = matrix_row_sub_index(game->payoffs[i][j], dbr_i);
        }
        matrix_copy_to_offset(C, A_ij, r_off, c_off);
        c_off += A_ij->ncols;
        //matrix_print(A_ij);
        //tC = matrix_augment(C, A_ij);
        //matrix_free(C);
        matrix_free(A_ij);
    }
}

matrix_t *construct_C1(polymatrix_t *game, strategy_profile_t *x, int *scount, int *K_x, int sKx, double delta)
{
    matrix_t **dbr = malloc(sizeof(matrix_t *) * sKx);
    int i, j, lp_rows = 0;
    for (i = 0; i < sKx; ++i) {
        dbr[i] = delta_best_response_i(game, x, delta, K_x[i]);
        lp_rows += dbr[i]->nrows;
    }

    if (lp_cols == 0) {
        lp_cols = sKx + 1;
        for (j = 0; j < game->players; ++j)
            lp_cols += scount[j];
    }

    matrix_t *C1 = matrix_alloc(lp_rows, lp_cols);

    int r_off = 0;
    for (i = 0; i < sKx; ++i) {
        //matrix_t *tC = construct_C1_i(game, scount, K_x[i], sKx, i, dbr[i], C1, r_off);
        construct_C1_i(game, scount, K_x[i], sKx, i, dbr[i], C1, r_off);
        //matrix_copy_to_offset(C1, tC, r_off, 0);
        r_off += dbr[i]->nrows;
        //matrix_free(tC);
        matrix_free(dbr[i]);
    }
    free(dbr);

    /*
    for (i = 1; i < sKx; ++i) {
        matrix_t *dbr_i = delta_best_response_i(game, x, delta, i);
        matrix_t *tC = construct_C1_i(game, scount, K_x[i], sKx, i, dbr_i);
        matrix_t *t =  matrix_augment_row(C1, tC);
        matrix_free(tC);
        matrix_free(C1);
        matrix_free(dbr_i);
        C1 = t;
    }*/
    return C1;
}

void construct_C2_i(polymatrix_t *game, strategy_profile_t *x, int *scount, int i, int sKx, int k, double *b, matrix_t *C2)
{
    int j;
    C2->data[k][0] = -1;
    C2->data[k][k + 1] = 1;
    matrix_t *vi = polymatrix_compute_vi(game, x, i);
    matrix_t *xt = matrix_trans(x->strategies[i]);
    matrix_t *u = matrix_mul_vec_mat(xt, vi);
    *b = -u->data[0][0];
    matrix_free(u);

    int c_off = sKx + 1;
    for (j = 0; j < game->players; ++j) {
        matrix_t *Ax;
        //printf("+++x++\n");
        //matrix_print(matrix_trans(x->strategies[j]));
        if (i == j) {
            Ax = matrix_trans(vi);
            matrix_mul_const_in(Ax, -1);
        } else if (game->payoffs[i][j]) {
            Ax = matrix_mul_vec_mat(xt, game->payoffs[i][j]);
            matrix_mul_const_in(Ax, -1);
        } else {
            Ax = matrix_alloc(1, scount[j]);
        }
        matrix_copy_to_offset(C2, Ax, i, c_off);
        c_off += Ax->ncols;
        matrix_free(Ax);
    }
    matrix_free(vi);
    matrix_free(xt);
}

matrix_t *construct_C2(polymatrix_t *game, strategy_profile_t *x, int *scount, int *K_x, int sKx, double *b)
{
    matrix_t *C2 = matrix_alloc(sKx, lp_cols);

    int i;
    for (i = 0; i < sKx; ++i)
        construct_C2_i(game, x, scount, K_x[i], sKx, i, &b[i], C2);

    return C2;
}

matrix_t *construct_C3(int sKx, int n, int *scount)
{
    int i, j, k;

    matrix_t *C3 = matrix_alloc(n, lp_cols);
    j = 1 + sKx;

    for (i = 0; i < n; ++i) {
        for(k = 0; k < scount[i]; ++k)
            C3->data[i][j++] = 1;
    }
    return C3;
}

cplp_t *form_descent_lp(polymatrix_t *game, strategy_profile_t *x, double delta, int *scount, int *K_x, int sKx)
{
    double *r = calloc(sKx, sizeof(double));
    matrix_t *C1 = construct_C1(game, x, scount, K_x, sKx, delta);
    matrix_t *C2 = construct_C2(game, x, scount, K_x, sKx, r);
    matrix_t *C3 = construct_C3(sKx, game->players, scount);
    matrix_t *A = matrix_augment_rows(3, C1, C2, C3);
    matrix_t *b = matrix_alloc(A->nrows, 1);
    int i;
    for (i = C1->nrows; i < b->nrows; ++i) {
        if (i - C1->nrows < C2->nrows)
            b->data[i][0] = r[i - C1->nrows];
        else
            b->data[i][0] = 1;
    }
    //printf("===============\n");
    //matrix_print(A);
    //matrix_print(b);
    matrix_t *c = matrix_alloc(A->ncols, 1);
    c->data[0][0] = 1;

    //matrix_print(c);
    cplp_t *lp = matrix_to_lp(A, b, c, 'L');
    lp->obj_sense = CPX_MIN;
    for (i = 0; i < C3->nrows; ++i)
        lp->sense[i + C1->nrows + C2->nrows] = 'E';

    lp->lb[0] = -CPX_INFBOUND;
    matrix_free(A);
    matrix_free(b);
    matrix_free(c);
    matrix_free(C1);
    matrix_free(C2);
    matrix_free(C3);
    free(r);

    return lp;
}

double delta_Df_i(polymatrix_t *game, strategy_profile_t *x, strategy_profile_t *xp, double delta, int i)
{
    matrix_t *v_i = polymatrix_compute_vi(game, x, i);
    matrix_t *vp_i = polymatrix_compute_vi(game, xp, i);

    //printf("v_i %d\n", i);
    //matrix_print(v_i);
    //printf("vp_i %d\n", i);
    //matrix_print(vp_i);
    //double Df_i = best_response_vi(vp_i, NULL);
    matrix_t *br_i = delta_best_response_vi(v_i, delta);
    //printf("===vi===\n");
    //matrix_print(v_i);
    //printf("====br==\n");
    //matrix_print(br_i);
    //printf("====vpi==\n");
    //matrix_print(vp_i);
    //printf("====xi===\n");
    //matrix_print(x->strategies[i]);
    //printf("====xpi===\n");
    //matrix_print(xp->strategies[i]);
    double Df_i = 0;

    int j;
    for (j = 0; j < br_i->nrows; ++j)
        Df_i = fmax(Df_i, vp_i->data[(int)br_i->data[j][0]][0]);

    //printf("Df_i %d = %lf ", i, Df_i);
    Df_i -= polymatrix_compute_ui_vi(vp_i, x, i);
    Df_i -= polymatrix_compute_ui_vi(v_i, xp, i);
    Df_i += polymatrix_compute_ui_vi(v_i, x, i);

    double t1 = polymatrix_compute_ui_vi(vp_i, x, i);
    double t2 = polymatrix_compute_ui_vi(v_i, xp, i);
    double t3 = polymatrix_compute_ui_vi(v_i, x, i);
    printf("- %lf - %lf + %lf", t1, t2, t3);

    return Df_i;
}

double delta_Df(polymatrix_t *game, strategy_profile_t *x, strategy_profile_t *xp, double delta, double eps, int *K_x, int sKx)
{
    //printf("================\n");
    double Df = -2;
    int i;
    for (i = 0; i < sKx; ++i){
        //printf("i %d\n", K_x[i]);
        double Df_i = delta_Df_i(game, x, xp, delta, K_x[i]);
        //printf(" %lf\n", Df_i);
        if (Df_i > Df)
            Df = Df_i;
    }

    //printf("Df = %lf - %lf = %lf\n",Df, eps, Df - eps);
    return Df - eps;
}

strategy_profile_t *solve_descent_lp(polymatrix_t *game, strategy_profile_t *x, int *scount, double delta, int *K_x, int m, double *Df)
{
    cplp_t *lp = form_descent_lp(game, x, delta, scount, K_x, m);
    cplp_sol_t *sol = cplp_solve(lp);
    *Df = sol->obj_val;
    //cplp_sol_print(sol);
    cplp_free(lp);

    strategy_profile_t *xp = strategy_alloc(game->players, scount);
    int i, j, k;
    k = m + 1;
    for (i = 0; i < game->players; ++i)
        for (j = 0; j < scount[i]; ++j)
            xp->strategies[i]->data[j][0] = sol->res[0][k++];

    //for (i = 0; i < game->players; ++i) {
        //matrix_print(matrix_trans(xp->strategies[i]));
        //matrix_print(matrix_trans(x->strategies[i]));
    //}

    cplp_sol_free(sol);
    return xp;
}

strategy_profile_t *new_point(strategy_profile_t *x, strategy_profile_t *xp, double eps)
{
    //printf("=======\n");
    strategy_profile_t *np = malloc(sizeof(strategy_profile_t));
    np->players = xp->players;
    np->strategies = malloc(sizeof(matrix_t *) * np->players);
    int i;
    for (i = 0; i < x->players; ++i) {
        matrix_t *tmp = matrix_sub(xp->strategies[i], x->strategies[i]);
        matrix_mul_const_in(tmp, eps);
        np->strategies[i] = matrix_add(x->strategies[i], tmp);
        matrix_free(tmp);
        //matrix_print(matrix_trans(xp->strategies[i]));
        //matrix_print(matrix_trans(x->strategies[i]));
        //matrix_print(matrix_trans(np->strategies[i]));
    }
    return np;
}

double find_eps(polymatrix_t *game, strategy_profile_t *x, strategy_profile_t *xp, double deps, int *scount)
{
    matrix_t *piece = piece_create(game, x, xp, scount);
    //fprintf(stderr,"%d %d\n", piece->nrows, piece->ncols);
    //print_piece(piece, stderr);
    //fprintf(stderr, "\n");
    //print_piece(piece_p(piece), stderr);
    struct timespec start, end, t;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);

    double steps = 1.0 / line_search;
    double eps, val, meps, mval;
    meps = deps;
    //mval = piece_eval(piece, 0);
    //printf("V0 -> %lf\n", mval);
    mval = piece_eval(piece, meps);
    //printf("Vt -> %lf\n", mval);
    double res, me, e;
    //strategy_profile_t *np;
    for (eps = 0; eps < 1; eps+=steps) {
        //np = new_point(x, xp, eps);
        //e = compute_eps(game, np, NULL);
        val = piece_eval(piece, eps);
        //printf("%lf %lf %lf\n", eps, val, e);
        if (val < mval) {
            mval = val;
            meps = eps;
        }
        /*
        if (e < me) {
            me = e;
            res = eps;
        }
        */
    }
    eps = 1;
    //np = new_point(x, xp, eps);
    val = piece_eval(piece, eps);
    if (val < mval) {
        mval = val;
        meps = eps;
    }
    //printf("%lf %lf %lf\n", eps, val, e);
    //exit(1);

    /*
    double me, res = 0;
    me = 1;
    //double eps;
    for (eps = 0; eps < 1; eps += 0.005) {
        strategy_profile_t *np = new_point(x, xp, eps);
        double e = compute_eps(game, np, NULL);
        if (e < me) {
            me = e;
            res = eps;
        }
        //fprintf(stderr, "%lf %lf\n", eps, e);
    }
    */

    //printf("Vm -> %lf %lf\n", meps, mval);
    //printf("Vm2 -> %lf  %lf\n", res, me);
    //return res;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end);
    double elapsed = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) * 1e-9;
    piece_time += elapsed;
    matrix_free(piece);
    return meps;
}

strategy_profile_t *descent(polymatrix_t *game, strategy_profile_t *xs, int *scount, double delta)
{
    double *eps_i = calloc(game->players, sizeof(double));
    double eps = compute_eps(game, xs, eps_i);
    g_eps = eps;
    int *K_x = calloc(game->players, sizeof(double));
    strategy_profile_t *x = xs;
    gx = x;
    double Df;
    double e1;
    double e = delta / (delta + 2);
    pivcount = 0;
    int m = get_K_x(eps, eps_i, game->players, K_x);

    //while(eps > 0.5 + delta) {
    while(1) {
        ++pivcount;
        //printf("-------------------------\n");
        //int m = get_K_x(eps, eps_i, game->players, K_x);
        //printf("SKx %d\n", m);
        strategy_profile_t *xp = solve_descent_lp(game, x, scount, delta, K_x, m, &Df);
        Df -= eps;

        strategy_profile_t *np;
        if (line_search > 0) {
            double dist = find_eps(game, x, xp, e, scount);
            np = new_point(x, xp, dist);
        } else {
            np = new_point(x, xp, e);
        }
        //printf("x\n");
        //matrix_print(matrix_trans(x->strategies[0]));
        //matrix_print(matrix_trans(x->strategies[1]));
        //printf("xp\n");
        //matrix_print(matrix_trans(xp->strategies[0]));
        //matrix_print(matrix_trans(xp->strategies[1]));
        //printf("np\n");
        //matrix_print(matrix_trans(np->strategies[0]));
        //matrix_print(matrix_trans(np->strategies[1]));
        //printf("Df %lf ", Df);
        //Df = delta_Df(game, x, xp, delta, eps, K_x, m);
        //printf("Eps %lf %lf\n", eps, Df);
        e1 = compute_eps(game, np, eps_i);
        //printf("Np %lf %lf\n", eps, e1);
        //e = delta / (delta + 2);
        if (err_le(eps, e1) && !err_eq(e1, eps)) {
            printf("Exiting %lf %lf\n", eps, e1);
            exit(1);
        }
        eps = e1;
        g_eps = eps;
        if (x != xs)
            strategy_free(x);
        strategy_free(xp);
        x = np;
        gx = x;
        if (Df >= -delta)
            break;
    }
    free(eps_i);
    free(K_x);
    return x;
}

strategy_profile_t *run_descent(polymatrix_t *game, double delta)
{
    int *scount = malloc(game->players * sizeof(int));
    count_strategies(game, scount);
    //polymatrix_normalize(game, scount);
    strategy_profile_t *x = strategy_alloc(game->players, scount);
    strategy_profile_t *xp  = descent(game, x, scount, delta);
    //double e1 = compute_eps(game, x, NULL);
    strategy_free(x);
    //strategy_print(xp);
    free(scount);
    return xp;
}

pthread_mutex_t lock = PTHREAD_MUTEX_INITIALIZER;
pthread_cond_t done = PTHREAD_COND_INITIALIZER;

struct descent_args {
    polymatrix_t *game;
    double delta;
};

void *solve_descent_thread(void *arg)
{
    int old_type;

    struct descent_args *args = (struct descent_args *)arg;
    pthread_setcanceltype(PTHREAD_CANCEL_ASYNCHRONOUS, &old_type);
    pthread_setcancelstate(PTHREAD_CANCEL_ENABLE, NULL);
    strategy_profile_t *x = run_descent(args->game, args->delta);

    pthread_cond_signal(&done);
    pthread_exit((void *)x);
}

void solve_descent(polymatrix_t *game, double delta)
{
    struct timespec start, end, t;

    clock_gettime(CLOCK_REALTIME, &t);
    t.tv_sec += TIMEOUT;
    pthread_t tid;
    
    struct descent_args args = {game, delta};
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);

    pthread_mutex_lock(&lock);
    pthread_create(&tid, NULL, solve_descent_thread, &args);
    int err = pthread_cond_timedwait(&done, &lock, &t);
    pthread_mutex_unlock(&lock);

    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end);
    double elapsed = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) * 1e-9;

    printf("Compute Time %lf\n", elapsed);
    printf("Iterations %d\n", pivcount);
    printf("Piece Time %lf\n", piece_time / pivcount);

    if (err != 0) {
        printf("Timeout\n");
        pthread_cancel(tid);
        printf("Eps %.20e\n", g_eps);
        strategy_print(gx);
        pthread_join(tid, NULL);
        return;
    }
    void *ptr;
    pthread_join(tid, &ptr);
    strategy_profile_t *x = (strategy_profile_t *)ptr;
    double e2 = compute_eps(game, x, NULL);
    printf("Eps %.20e\n", e2);
    strategy_print(x);
    strategy_free(x);
}
