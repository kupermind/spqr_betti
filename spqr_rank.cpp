#include <stdio.h>
#include <vector>
#include <set>
#include <utility>
#include <map>
#include <fstream>

#include "spqr.hpp"
#include "SuiteSparseQR.hpp"

int main(int argc, char * argv[])
{
    Long m, n, i, j, li = 1;
    double val;
    cholmod_sparse *A = new cholmod_sparse;
    int let_transpose = 0;
    if (argc == 2)
        let_transpose = atoi(argv[1]);


    /// Initial set up
    A->packed = TRUE;
    A->sorted = FALSE;
    A->nz = NULL;
    A->itype = CHOLMOD_LONG;
    A->dtype = CHOLMOD_DOUBLE;
    A->stype = 0;
    A->z = NULL;
    A->xtype = CHOLMOD_REAL;

    /// Getting the number of rows, columns
    std::ifstream f(argv[1]);
    f >> m >> n;

    /// The matrix is column-driven
    if (0)
//    if (m <= n || !let_transpose)
    {
        A->nrow = m;
        A->ncol = n;
        printf("Loading %ld by %ld matrix and converting to column-driven format...\n", A->nrow, A->ncol);
    }
    else
    {
        A->nrow = n;
        A->ncol = m;
        printf("Loading %ld by %ld matrix...\n", A->nrow, A->ncol);
    }
    fflush(stdout);

    std::vector<Long> irarr;
    std::vector<Long> jcarr(A->ncol + 1, 0);
    std::vector<double> prarr;

    /// Getting the rest of arrays
    /// If m > n, then do the transpose using the row information
    if (1)
//    if (m >= n && let_transpose)
    {
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
//            jcarr.push_back(j);
            prarr.push_back(val);        
        }
    }
    else
    /// Otherwise, convert into the column-driven format
    {
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
    }
    f.close();
    irarr.resize(irarr.size());
    prarr.resize(prarr.size());
   
    printf("Loaded the matrix, calculating rank...\n");
    fflush(stdout);

    A->p = &jcarr[0];
    A->i = &irarr[0];
    A->x = &prarr[0];
    Long *Ap;
    Ap = (Long*)A->p;
    A->nzmax = Ap[A->ncol];

    /// Analyze, factorize, and get the results
    int order = 6;
//    int order = 0;
    double tol = -2;
    Long econ = A->nrow;
    cholmod_sparse *Q, *R;
    Long *E;
    cholmod_common *cc = new cholmod_common;
    cholmod_l_defaults(cc);
    cholmod_l_start(cc);

//    Long rank = SuiteSparseQR <double> (order, tol, econ, A, &Q, &R, &E, cc);
    Long rank = SuiteSparseQR <double> (order, tol, econ, A, &R, &E, cc);

    printf("The rank is %ld\n", rank);
    

    delete cc;
    delete A;
}
