#include <stdio.h>
#include <vector>
#include <set>
#include <utility>
#include <map>
#include <fstream>

#include "spqr.hpp"
#include "SuiteSparseQR.hpp"


void mat_setup(cholmod_sparse *A)
{
    /// Initial set up
    A->packed = TRUE;
    A->sorted = FALSE;
    A->nz = NULL;
    A->itype = CHOLMOD_LONG;
    A->dtype = CHOLMOD_DOUBLE;
    A->stype = 0;
    A->z = NULL;
    A->xtype = CHOLMOD_REAL;
}


void load_normal(cholmod_sparse *A, char* filename)
{
    Long m, n, i, j;
    double val;

    /// Getting the number of rows, columns
    std::ifstream f(filename);
    f >> m >> n;

    A->nrow = m;
    A->ncol = n;
    printf("Loading %ld by %ld matrix and converting to column-driven format...\n", A->nrow, A->ncol);

    std::vector<Long> irarr;
    std::vector<Long> jcarr(A->ncol + 1, 0);
    std::vector<double> prarr;
    /// Convert into the column-driven format
    A->sorted = TRUE;

    /// First, build a map of the matrix
    std::map<Long, std::set<std::pair<Long, double> > > mat_map;
    while (!f.eof())
    {
        f >> i >> j >> val;
        --i; --j;
        Long jcnt = mat_map.count(j);
        if (!jcnt)
        {
            std::set<std::pair<Long, double> > init_set;
            mat_map[j] = init_set;
        }
        mat_map[j].insert(std::make_pair(i, val));
    }

    /// Then construct all the needed arrays
    std::map<Long, std::set<std::pair<Long, double> > >::iterator map_it = mat_map.begin();
    for (i=1; map_it != mat_map.end(); ++map_it, ++i)
    {
        jcarr[i] = map_it->second.size() + jcarr[i-1];
        std::set<std::pair<Long, double> >::iterator set_it;
        set_it = map_it->second.begin();
        for (; set_it != map_it->second.end(); ++set_it)
        {
            irarr.push_back(set_it->first);
            prarr.push_back(set_it->second);
        }
    }

    mat_map.clear();
    f.close();
    irarr.resize(irarr.size());
    prarr.resize(prarr.size());

    Long *ap = new Long[jcarr.size()];
    A->p = (void*)ap;
    std::copy(jcarr.begin(), jcarr.end(), ap);
    Long *ai = new Long[irarr.size()];
    A->i = (void*)ai;
    std::copy(irarr.begin(), irarr.end(), ai);
    double *ax = new double[prarr.size()];
    A->x = (void*)ax;
    std::copy(prarr.begin(), prarr.end(), ax);
    Long *Ap;
    Ap = (Long*)A->p;
    A->nzmax = Ap[A->ncol];
}

void load_transpose(cholmod_sparse *A, char* filename)
{
    Long m, n, i, j, li = 1;
    double val;

    /// Getting the number of rows, columns
    std::ifstream f(filename);
    f >> m >> n;

    A->nrow = n;
    A->ncol = m;
    printf("Loading transposed %ld by %ld matrix...\n", A->nrow, A->ncol);

    std::vector<Long> irarr;
    std::vector<Long> jcarr(A->ncol + 1, 0);
    std::vector<double> prarr;

    /// Getting the rest of arrays
    while (!f.eof())
    {
        f >> i >> j >> val;
        i; --j;
        if(li != i)
        {
            li = i;
            jcarr[i] = jcarr[i-1];
        }
        ++jcarr[i];

        irarr.push_back(j);
        prarr.push_back(val);
    }

    f.close();
    irarr.resize(irarr.size());
    prarr.resize(prarr.size());

    Long *ap = new Long[jcarr.size()];
    A->p = (void*)ap;
    std::copy(jcarr.begin(), jcarr.end(), ap);
    Long *ai = new Long[irarr.size()];
    A->i = (void*)ai;
    std::copy(irarr.begin(), irarr.end(), ai);
    double *ax = new double[prarr.size()];
    A->x = (void*)ax;
    std::copy(prarr.begin(), prarr.end(), ax);
    Long *Ap;
    Ap = (Long*)A->p;
    A->nzmax = Ap[A->ncol];
}


int main(int argc, char * argv[])
{
    cholmod_sparse *AT = new cholmod_sparse;
    cholmod_sparse *BT = new cholmod_sparse;
    cholmod_common *cc = new cholmod_common;
    cholmod_l_defaults(cc);
    cholmod_l_start(cc);

    /// Loading the matrices
    mat_setup(AT);
    mat_setup(BT);

    load_transpose(AT, argv[1]);
    cholmod_sort(AT, cc);
//    load_normal(A, argv[1]);
    load_transpose(BT, argv[2]);
    cholmod_sort(BT, cc);

    printf("Loaded the matrices\n");
    fflush(stdout);

    /// Calculating Laplacian
    cholmod_sparse *A = cholmod_l_transpose(AT, 1, cc);
//printf("A: nrow %ld, ncol %ld, nzmax %ld, packed %d, sorted %d, stype %d, xtype %d, first elem %lf\n", A->nrow, A->ncol, A->nzmax, A->packed, A->sorted, A->stype, A->xtype, val_arr[0]);
//FILE *fa = fopen("sparseA.txt", "w");
//cholmod_l_write_sparse(fa, A, NULL, NULL, cc);
//fclose(fa);
//printf("Transposed the matrix A, with dimensions %ld by %ld\n", A->nrow, A->ncol);
//fflush(stdout);

    cholmod_sparse *ATA = cholmod_l_ssmult(AT, A, 1, 1, 1, cc);
//printf("Multiplied AT * A: %ld by %ld\n", ATA->nrow, ATA->ncol);
//fflush(stdout);

    cholmod_sparse *B = cholmod_l_transpose(BT, 1, cc);
//printf("Transposed the matrix B, with dimensions %ld by %ld\n", B->nrow, B->ncol);
//fflush(stdout);

//    cholmod_sparse *BBT = cholmod_l_ssmult(B, BT, 0, 1, 1, cc);
    cholmod_sparse *BBT = cholmod_l_ssmult(B, BT, 1, 1, 1, cc);
//printf("Multiplied B * BT: %ld by %ld\n", BBT->nrow, BBT->ncol);
//fflush(stdout);


//    cholmod_sparse *ATA = cholmod_l_aat(A, NULL, 0, 1, cc);
    delete A, AT, B, BT;

    printf("Transposed for building the Laplacian\n");
    fflush(stdout);

    if (ATA->nrow != BBT->nrow)
    {
        printf("ERROR: sizes of A'*A and B*B' are different! Aborting...\n");
        abort();
    }
    double alpha[2] = {1, 1};
    double beta[2] = {1, 1};
    cholmod_sparse *L = cholmod_l_add(ATA, BBT, alpha, beta, 1, 1, cc);

    delete ATA, BBT;

    printf("Built the Laplacian %ld by %ld, calculating Betti number...\n", L->nrow, L->ncol);
    fflush(stdout);

    /// Analyze, factorize, and get the results
    int order = 6;
//    int order = 4;
    double tol = -2;
    Long econ = L->nrow;
    cholmod_sparse *Q, *R;
    Long *E;

    Long rank = SuiteSparseQR <double> (order, tol, econ, L, &R, &E, cc);

    printf("The Betti number is %ld\n", L->nrow - rank);
    

    delete cc;
    delete L;

}
