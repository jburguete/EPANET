/**
 * \file genmmd.h
 * \brief Header file for multiple minimum degree row re-ordering algorithm.
 * \authors see AUTHORS.
 * \copyright see AUTHORS.
 */


#ifndef GENMMD__H
#define GENMMD__H 1

void
genmmd (int *neqns, int *xadj, int *adjncy, int *invp, int *perm,
        int *delta, int *dhead, int *qsize, int *llist, int *marker,
        int *maxint, int *nofsub);

#endif
