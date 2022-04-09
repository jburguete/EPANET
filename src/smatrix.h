/**
 * \file smatrix.h
 * \brief Header file to solve a sparse set of linear equations.
 * \authors see AUTHORS.
 * \copyright see AUTHORS.
 * 
 * This file contains the sparse matrix routines used to solve a network's
 * hydraulic equations. The functions exported by this module are:
 * createsparse() -- called from openhyd() in hydraul.c,
 * freesparse()   -- called from closehyd() in hydraul.c,
 * linsolve()     -- called from hydsolve() in hydsolver.c.
 * 
*/

#ifndef SMATRIX__H
#define SMATRIX__H 1

int createsparse (Network * net, Smatrix * sm);

/**
 * function to free memory used for sparse matrix storage.
 */
static inline void
freesparse (Smatrix * sm)       ///< sparse matrix struct.
{

  FREE (sm->Order);
  FREE (sm->Row);
  FREE (sm->Ndx);
  FREE (sm->XLNZ);
  FREE (sm->NZSUB);
  FREE (sm->LNZ);
  FREE (sm->Aij);
  FREE (sm->Aii);
  FREE (sm->F);
  FREE (sm->temp);
  FREE (sm->link);
  FREE (sm->first);
}

/**
 * function to solve sparse symmetric system of linear equations using Cholesky
 * factorization.
 * Output:  sm->F = solution values
 *
 * \return 0 if solution found, or index of equation causing system to be
 * ill-conditioned.
 *
 * NOTE: This procedure assumes that the solution matrix has been
 * symbolically factorized with the positions of the lower triangular,
 * off-diagonal, non-zero coeffs. stored in the following integer arrays:
 * 
 * XLNZ  (start position of each column in NZSUB),
 * 
 * NZSUB (row index of each non-zero in each column),
 * 
 * LNZ   (position of each NZSUB entry in Aij array),
 *
 * This procedure has been adapted from subroutines GSFCT and GSSLV in the book
 * "Computer Solution of Large Sparse Positive Definite Systems" by A. George
 * and J. W-H Liu (Prentice-Hall, 1981).
*/
static inline int
linsolve (Smatrix * sm,         ///< sparse matrix struct.
          int n)                ///< number of equations.
{
  double *Aii = sm->Aii;
  double *Aij = sm->Aij;
  double *B = sm->F;
  double *temp = sm->temp;
  int *LNZ = sm->LNZ;
  int *XLNZ = sm->XLNZ;
  int *NZSUB = sm->NZSUB;
  int *link = sm->link;
  int *first = sm->first;

  int i, istop, istrt, isub, j, k, kfirst, newk;
  double bj, diagj, ljk;

  memset (temp, 0, (n + 1) * sizeof (double));
  memset (link, 0, (n + 1) * sizeof (int));
  memset (first, 0, (n + 1) * sizeof (int));

  // Begin numerical factorization of matrix A into L
  //   Compute column L(*,j) for j = 1,...n
  for (j = 1; j <= n; j++)
    {
      // For each column L(*,k) that affects L(*,j):
      diagj = 0.0;
      newk = link[j];
      k = newk;
      while (k != 0)
        {
          // Outer product modification of L(*,j) by
          // L(*,k) starting at first[k] of L(*,k)
          newk = link[k];
          kfirst = first[k];
          ljk = Aij[LNZ[kfirst]];
          diagj += ljk * ljk;
          istrt = kfirst + 1;
          istop = XLNZ[k + 1] - 1;
          if (istop >= istrt)
            {

              // Before modification, update vectors 'first'
              // and 'link' for future modification steps
              first[k] = istrt;
              isub = NZSUB[istrt];
              link[k] = link[isub];
              link[isub] = k;

              // The actual mod is saved in vector 'temp'
              for (i = istrt; i <= istop; i++)
                {
                  isub = NZSUB[i];
                  temp[isub] += Aij[LNZ[i]] * ljk;
                }
            }
          k = newk;
        }

      // Apply the modifications accumulated
      // in 'temp' to column L(*,j)
      diagj = Aii[j] - diagj;
      if (diagj <= 0.0)         // Check for ill-conditioning
        {
          return j;
        }
      diagj = sqrt (diagj);
      Aii[j] = diagj;
      istrt = XLNZ[j];
      istop = XLNZ[j + 1] - 1;
      if (istop >= istrt)
        {
          first[j] = istrt;
          isub = NZSUB[istrt];
          link[j] = link[isub];
          link[isub] = j;
          for (i = istrt; i <= istop; i++)
            {
              isub = NZSUB[i];
              bj = (Aij[LNZ[i]] - temp[isub]) / diagj;
              Aij[LNZ[i]] = bj;
              temp[isub] = 0.0;
            }
        }
    }                           // next j

  // Foward substitution
  for (j = 1; j <= n; j++)
    {
      bj = B[j] / Aii[j];
      B[j] = bj;
      istrt = XLNZ[j];
      istop = XLNZ[j + 1] - 1;
      if (istop >= istrt)
        {
          for (i = istrt; i <= istop; i++)
            {
              isub = NZSUB[i];
              B[isub] -= Aij[LNZ[i]] * bj;
            }
        }
    }

  // Backward substitution
  for (j = n; j >= 1; j--)
    {
      bj = B[j];
      istrt = XLNZ[j];
      istop = XLNZ[j + 1] - 1;
      if (istop >= istrt)
        {
          for (i = istrt; i <= istop; i++)
            {
              isub = NZSUB[i];
              bj -= Aij[LNZ[i]] * B[isub];
            }
        }
      B[j] = bj / Aii[j];
    }
  return 0;
}

#endif
