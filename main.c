#include <stdio.h>
#include <stdlib.h>
#include "../util/polymatrix.h"
#include "../util/io.h"
#include "descent.h"
#include <time.h>

int main(int argc, char **argv)
{
    FILE *f = fopen(argv[1], "r");
    polymatrix_t *game = read_polymatrix_from_file(f, NULL,NULL);
    //polymatrix_print(game);
    int *scount = malloc(game->players * sizeof(int));
    count_strategies(game, scount);

    //polymatrix_print(game);
    struct timespec start, end;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);
    polymatrix_normalize(game, scount);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end);
    double elapsed = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) * 1e-9;
    printf("Normal Time %lf\n", elapsed);
    //polymatrix_print(game);
    //clock_t start = clock(), diff;
    if (argc > 3)
        line_search = atof(argv[3]);
    solve_descent(game, atof(argv[2]));
    //diff = clock() - start;
    //printf("%d\n", (int)(diff / CLOCKS_PER_SEC));
    printf("\n");
    polymatrix_free(game);
    free(scount);
    return 0;
}
