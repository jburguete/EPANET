/**
 * \file smatrix.c
 * \brief Source file to solve a sparse set of linear equations.
 * \authors see AUTHORS.
 * \copyright see AUTHORS.
 * 
 * This file contains the sparse matrix routines used to solve a network's
 * hydraulic equations. The functions exported by this module are:
 *
 * createsparse() -- called from openhyd() in hydraul.c,
 * 
 * freesparse()   -- called from closehyd() in hydraul.c,
 * 
 * linsolve()     -- called from hydsolve() in hydsolver.c.
 * 
*/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <time.h>  //For optional timer macros
#include <glib.h>

#include "text.h"
#include "types.h"
#include "funcs.h"
#include "smatrix.h"

// The multiple minimum degree re-ordering routine (see genmmd.c)
extern int genmmd(int *neqns, int *xadj, int *adjncy, int *invp, int *perm,
                  int *delta, int *dhead, int *qsize, int *llist, int *marker,
                  int *maxint, int *nofsub);

/**
 * function to allocate memory for representing a sparse matrix.
 *
 * \return error code.
 */
static inline int  allocsmatrix(Smatrix *sm,   ///< sparse matrix struct.
                                int Nnodes, int Nlinks)
{
    int errcode = 0;

    // Memory for linear eqn. solver allocated in alloclinsolve().
    sm->Aij   = NULL;
    sm->Aii   = NULL;
    sm->F     = NULL;
    sm->temp  = NULL;
    sm->link  = NULL;
    sm->first = NULL;

    // Memory for representing sparse matrix data structure
    sm->Order  = (int *) calloc(Nnodes+1,  sizeof(int));
    sm->Row    = (int *) calloc(Nnodes+1,  sizeof(int));
    sm->Ndx    = (int *) calloc(Nlinks+1,  sizeof(int));
    ERRCODE(MEMCHECK(sm->Order));
    ERRCODE(MEMCHECK(sm->Row));
    ERRCODE(MEMCHECK(sm->Ndx));
    return errcode;
}

/**
 * function to allocate memory used by linear eqn. solver.
 *
 * \return error code.
 */
static inline int  alloclinsolve(Smatrix *sm,   ///< sparse matrix struct.
                                 int n)
{
    int errcode = 0;
    n = n + 1;    // All arrays are 1-based

    sm->Aij   = (double *)calloc(sm->Ncoeffs + 1, sizeof(double));
    sm->Aii   = (double *)calloc(n, sizeof(double));
    sm->F     = (double *)calloc(n, sizeof(double));
    sm->temp  = (double *)calloc(n, sizeof(double));
    sm->link  = (int *)calloc(n, sizeof(int));
    sm->first = (int *)calloc(n, sizeof(int));
    ERRCODE(MEMCHECK(sm->Aij));
    ERRCODE(MEMCHECK(sm->Aii));
    ERRCODE(MEMCHECK(sm->F));
    ERRCODE(MEMCHECK(sm->temp));
    ERRCODE(MEMCHECK(sm->link));
    ERRCODE(MEMCHECK(sm->first));
    return errcode;
}

/**
 * function to remove parallel links from nodal adjacency lists.
 */
static inline void  xparalinks(Padjlist *Adjlist, int Nnodes)
{
    int    i;
    Padjlist    alink,       // Current item in adjacency list
                blink;       // Previous item in adjacency list

    // Scan adjacency list of each node
    for (i = 1; i <= Nnodes; i++)
    {
        alink = Adjlist[i];               // First item in list
        blink = NULL;
        while (alink != NULL)
        {
            if (alink->node == 0)              // Parallel link marker found
            {
                if (blink == NULL)             // This holds at start of list
                {
                    Adjlist[i] = alink->next;
                    free(alink);                // Remove item from list
                    alink = Adjlist[i];
                }
                else                           // This holds for interior of list
                {
                    blink->next = alink->next;
                    free(alink);                // Remove item from list
                    alink = blink->next;
                }
            }
            else
            {
                blink = alink;                // Move to next item in list
                alink = alink->next;
            }
        }
    }
}

/**
 * function to check for parallel links between nodes i and j.
 *
 * \return 1 if link k parallels another link, else 0.
 */
static inline int  paralink(Padjlist *Adjlist,
                            Smatrix *sm,   ///< sparse matrix struct.
                            int i,      ///< index of start node of link.
                            int j,      ///< index of end node of link.
                            int k)      ///< link index.
{
    Padjlist alink;
    for (alink = Adjlist[i]; alink != NULL; alink = alink->next)
    {
        // Link || to k (same end nodes)
        if (alink->node == j)
        {
            // Assign Ndx entry to this link
            sm->Ndx[k] = alink->link;
            return(1);
        }
    }
    // Ndx entry if link not parallel
    sm->Ndx[k] = k;
    return(0);
}

/**
 * function to build linked list of non-parallel links adjacent to each node.
 *
 * \return error code.
 */
static inline int  localadjlists(Network *net,
                                 Smatrix *sm)   ///< sparse matrix struct.
{
    int    i, j, k;
    int    pmark = 0;     // parallel link marker
    int    errcode = 0;
    Padjlist  alink;

    // Create an array of adjacency lists
    freeadjlists(net);
    net->Adjlist = (Padjlist *)calloc(net->Nnodes + 1, sizeof(Padjlist));
    if (net->Adjlist == NULL) return 101;

    // For each link, update adjacency lists of its end nodes
    for (k = 1; k <= net->Nlinks; k++)
    {
        i = net->Link[k].N1;
        j = net->Link[k].N2;
        pmark = paralink(net->Adjlist, sm, i, j, k);  // Parallel link check

        // Include link in start node i's list
        alink = (struct Sadjlist *) malloc(sizeof(struct Sadjlist));
        if (alink == NULL) return(101);
        if (!pmark) alink->node = j;
        else        alink->node = 0;         // Parallel link marker
        alink->link = k;
        alink->next = net->Adjlist[i];
        net->Adjlist[i] = alink;

        // Include link in end node j's list
        alink = (struct Sadjlist *) malloc(sizeof(struct Sadjlist));
        if (alink == NULL) return(101);
        if (!pmark) alink->node = i;
        else        alink->node = 0;         // Parallel link marker
        alink->link = k;
        alink->next = net->Adjlist[j];
        net->Adjlist[j] = alink;
    }

    // Remove parallel links from adjacency lists
    xparalinks(net->Adjlist, net->Nnodes);
    return errcode;
}

/**
 * function to re-order nodes to minimize # of non-zeros that will appear in
 * factorized solution matrix.
 *
 * \return 1 if successful, 0 if not.
 */
static inline int reordernodes(Network *net,
                               Smatrix *sm)   ///< sparse matrix struct.
{
    int k, knode, m, njuncs, nlinks;
    int delta = -1;
    int nofsub = 0;
    int maxint = INT_MAX;   //defined in limits.h
    int errcode;
    Padjlist alink;

    // Local versions of node adjacency lists
    int *adjncy = NULL;
    int *xadj   = NULL;

    // Work arrays
    int *dhead = NULL;
    int *qsize = NULL;
    int *llist = NULL;
    int *marker = NULL;

    // Default ordering
    for (k = 1; k <= net->Nnodes; k++)
    {
        sm->Row[k] = k;
        sm->Order[k] = k;
    }
    njuncs = net->Njuncs;
    nlinks = net->Nlinks;

    // Allocate memory
    adjncy = (int *) calloc(2*nlinks+1, sizeof(int));
    xadj   = (int *) calloc(njuncs+2, sizeof(int));
    dhead  = (int *) calloc(njuncs+1, sizeof(int));
    qsize  = (int *) calloc(njuncs + 1, sizeof(int));
    llist  = (int *) calloc(njuncs + 1, sizeof(int));
    marker = (int *) calloc(njuncs + 1, sizeof(int));
    if (adjncy && xadj && dhead && qsize && llist && marker)
    {
        // Create local versions of node adjacency lists
        xadj[1] = 1;
        m = 1;
        for (k = 1; k <= njuncs; k++)
        {
            for (alink = net->Adjlist[k]; alink != NULL; alink = alink->next)
            {
                knode = alink->node;
                if (knode > 0 && knode <= njuncs)
                {
                    adjncy[m] = knode;
                    m++;
                }
            }
            xadj[k+1] = m;
        }

        // Generate a multiple minimum degree node re-ordering
        genmmd(&njuncs, xadj, adjncy, sm->Row, sm->Order, &delta,
               dhead, qsize, llist, marker, &maxint, &nofsub);
        errcode = 0;
    }
    else errcode = 101;  //insufficient memory

    // Free memory
    FREE(adjncy);
    FREE(xadj);
    FREE(dhead);
    FREE(qsize);
    FREE(llist);
    FREE(marker);
    return errcode;
}

/**
 * function to check if nodes i and j are already linked.
 *
 * \return 1 if nodes i and j are linked, 0 if not.
*/
static inline int linked(Padjlist * Adjlist,
                         int i, ///< node index
                         int j) ///< node index
{
    Padjlist alink;
    for (alink = Adjlist[i]; alink != NULL; alink = alink->next)
    {
        if (alink->node == j) return 1;
    }
    return 0;
}

/**
 * function to augment node i's adjacency list with node j.
 *
 * \return 1 if successful, 0 if not.
*/
static int  addlink(Padjlist * Adjlist,
                    int i,      ///< node index.
                    int j,      ///< node index.
                    int n)      ///< link index.
{
    Padjlist alink;
    alink = (struct Sadjlist *) malloc(sizeof(struct Sadjlist));
    if (alink == NULL) return 0;
    alink->node = j;
    alink->link = n;
    alink->next = Adjlist[i];
    Adjlist[i] = alink;
    return 1;
}

/**
 * function to link end of current adjacent link to end nodes of all links that
 * follow it on adjacency list.
 *
 * \return 1 if successful, 0 if not.
 */
static inline int  newlink(Padjlist *Adjlist,
                           Smatrix *sm,   ///< sparse matrix struct.
                           Padjlist alink)
                           ///< element of node's adjacency list.
{
    int inode, jnode;
    Padjlist blink;

    // Scan all entries in adjacency list that follow anode.
    inode = alink->node;             // End node of connection to anode
    for (blink = alink->next; blink != NULL; blink = blink->next)
    {
        jnode = blink->node;          // End node of next connection

        // If jnode still active, and inode not connected to jnode,
        // then add a new connection between inode and jnode.
        if (jnode > 0 && sm->Degree[jnode] > 0)  // jnode still active
        {
            if (!linked(Adjlist, inode, jnode))      // inode not linked to jnode
            {
                // Since new connection represents a non-zero coeff.
                // in the solution matrix, update the coeff. count.
                sm->Ncoeffs++;

                // Update adjacency lists for inode & jnode to
                // reflect the new connection.
                if (!addlink(Adjlist, inode, jnode, sm->Ncoeffs)) return 0;
                if (!addlink(Adjlist, jnode, inode, sm->Ncoeffs)) return 0;
                sm->Degree[inode]++;
                sm->Degree[jnode]++;
            }
        }
    }
    return 1;
}

/**
 * function to create new entries in knode's adjacency list for all unlinked
 * pairs of active nodes that are adjacent to knode.
 *
 * \return 1 if successful, 0 if not.
 */
static inline int  growlist(Network * net,
                            Smatrix * sm,   ///< sparse matrix struct.
                            int knode)  ///< node index.
{

    int node;
    Padjlist alink;

    // Iterate through all nodes connected to knode
    for (alink = net->Adjlist[knode]; alink != NULL; alink = alink -> next)
    {
        node = alink->node;                   // End node of connecting link
        if (node > 0 && sm->Degree[node] > 0) // End node is active
        {
            sm->Degree[node]--;           // Reduce degree of adjacency
            if (!newlink(net->Adjlist, sm, alink))      // Add to adjacency list
            {
                return 0;
            }
        }
  }
  return 1;
}

/**
 * function to symbolically factorize the solution matrix in terms of its
 * adjacency lists.
 *
 * \return error code.
*/
static inline int factorize(Network * net,
                            Smatrix * sm)   ///< sparse matrix struct.
{
    int k, knode;
    int errcode = 0;
    Padjlist alink;

    // Find degree of each junction node
    sm->Degree = (int *)calloc(net->Nnodes + 1, sizeof(int));
    if (sm->Degree == NULL) return 101;

    // NOTE: For purposes of node re-ordering, Tanks (nodes with
    //       indexes above Njuncs) have zero degree of adjacency.

    for (k = 1; k <= net->Njuncs; k++)
    {
        for (alink = net->Adjlist[k]; alink != NULL; alink = alink->next)
        {
            if (alink->node > 0) sm->Degree[k]++;
        }
    }

    // Augment each junction's adjacency list to account for
    // new connections created when solution matrix is solved.
    // NOTE: Only junctions (indexes <= Njuncs) appear in solution matrix.
    for (k = 1; k <= net->Njuncs; k++)          // Examine each junction
    {
        knode = sm->Order[k];                   // Re-ordered index
        if (!growlist(net, sm, knode))               // Augment adjacency list
        {
            errcode = 101;
            break;
        }
        sm->Degree[knode] = 0;                  // In-activate node
    }
    free(sm->Degree);
    return errcode;
}

/**
 * function to store row indexes of non-zeros of each column of lower triangular
 * portion of factorized matrix.
 *
 * \return error code.
*/
static inline int  storesparse(Network * net,
                               Smatrix * sm,   ///< sparse matrix struct.
                               int n)   ///< number of rows in solution matrix.
{
    int i, ii, j, k, l, m;
    int errcode = 0;
    Padjlist alink;

    // Allocate sparse matrix storage
    sm->XLNZ  = (int *) calloc(n+2, sizeof(int));
    sm->NZSUB = (int *) calloc(sm->Ncoeffs+2, sizeof(int));
    sm->LNZ   = (int *) calloc(sm->Ncoeffs+2, sizeof(int));
    ERRCODE(MEMCHECK(sm->XLNZ));
    ERRCODE(MEMCHECK(sm->NZSUB));
    ERRCODE(MEMCHECK(sm->LNZ));
    if (errcode) return errcode;

    // Generate row index pointers for each column of matrix
    k = 0;
    sm->XLNZ[1] = 1;
    for (i = 1; i <= n; i++)            // column
    {
        m = 0;
        ii = sm->Order[i];
        for (alink = net->Adjlist[ii]; alink != NULL; alink = alink->next)
        {
            if (alink->node == 0) continue;
            j = sm->Row[alink->node];    // row
            l = alink->link;
            if (j > i && j <= n)
            {
                m++;
                k++;
                sm->NZSUB[k] = j;
                sm->LNZ[k] = l;
            }
        }
        sm->XLNZ[i+1] = sm->XLNZ[i] + m;
    }
    return errcode;
}

/**
 * function to determine sparse storage scheme for transpose of a matrix.
*/
static void transpose(int n,    ///< matrix order.
                      int *il,
                      int *jl,
                      int *xl,  ///< sparse storage scheme for original matrix.
                      int *ilt,
                      int *jlt,
                      int *xlt,
                      ///< sparse storage scheme for transposed matrix.
                      int *nzt)       ///< work array.
{
    int  i, j, k, kk;

    for (i = 1; i <= n; i++) nzt[i] = 0;
    for (i = 1; i <= n; i++)
    {
        for (k = il[i]; k < il[i+1]; k++)
        {
            j = jl[k];
            kk = ilt[j] + nzt[j];
            jlt[kk] = i;
            xlt[kk] = xl[k];
            nzt[j]++;
        }
    }
}

/**
 * function to put row indexes in ascending order in NZSUB.
 *
 * \return eror code.
 */
static inline int sortsparse(Smatrix *sm,   ///< sparse matrix struct.
                             int n)     ///< number of rows in solution matrix.
{
    int  i, k;
    int  *xlnzt, *nzsubt, *lnzt, *nzt;
    int  errcode = 0;

    int *LNZ = sm->LNZ;
    int *XLNZ = sm->XLNZ;
    int *NZSUB = sm->NZSUB;

    xlnzt  = (int *) calloc(n+2, sizeof(int));
    nzsubt = (int *) calloc(sm->Ncoeffs+2, sizeof(int));
    lnzt   = (int *) calloc(sm->Ncoeffs+2, sizeof(int));
    nzt    = (int *) calloc(n+2, sizeof(int));
    ERRCODE(MEMCHECK(xlnzt));
    ERRCODE(MEMCHECK(nzsubt));
    ERRCODE(MEMCHECK(lnzt));
    ERRCODE(MEMCHECK(nzt));
    if (!errcode)
    {
        // Count # non-zeros in each row
        for (i = 1; i <= n; i++) nzt[i] = 0;
        for (i = 1; i <= n; i++)
        {
            for (k = XLNZ[i]; k < XLNZ[i+1]; k++) nzt[NZSUB[k]]++;
        }
        xlnzt[1] = 1;
        for (i = 1; i <= n; i++) xlnzt[i+1] = xlnzt[i] + nzt[i];

        // Transpose matrix twice to order column indexes
        transpose(n, XLNZ, NZSUB, LNZ, xlnzt, nzsubt, lnzt, nzt);
        transpose(n, xlnzt, nzsubt, lnzt, XLNZ, NZSUB, LNZ, nzt);
    }

    // Reclaim memory
    free(xlnzt);
    free(nzsubt);
    free(lnzt);
    free(nzt);
    return errcode;
}

/**
 * function to create sparse representation of coeff. matrix.
 *
 * \return error code.
 */
int  createsparse(Network * net,
                  Smatrix * sm)   ///< sparse matrix struct.
{
    int errcode = 0;

//    cleartimer(SmatrixTimer);
//    starttimer(SmatrixTimer);

    // Allocate sparse matrix data structures
    errcode = allocsmatrix(sm, net->Nnodes, net->Nlinks);
    if (errcode) return errcode;

    // Build a local version of node-link adjacency lists
    // with parallel links removed
    errcode = localadjlists(net, sm);
    if (errcode) return errcode;

    // Re-order nodes to minimize number of non-zero coeffs.
    // in factorized solution matrix
    ERRCODE(reordernodes(net, sm));

    // Factorize solution matrix by updating adjacency lists
    // with non-zero connections due to fill-ins
    sm->Ncoeffs = net->Nlinks;
    ERRCODE(factorize(net, sm));

    // Allocate memory for sparse storage of positions of non-zero
    // coeffs. and store these positions in vector NZSUB
    ERRCODE(storesparse(net, sm, net->Njuncs));

    // Free memory used for local adjacency lists and sort
    // row indexes in NZSUB to optimize linsolve()
    freeadjlists(net);
    ERRCODE(sortsparse(sm, net->Njuncs));

    // Allocate memory used by linear eqn. solver
    ERRCODE(alloclinsolve(sm, net->Nnodes));

    // Re-build adjacency lists for future use
    ERRCODE(buildadjlists(net));
    return errcode;
}
