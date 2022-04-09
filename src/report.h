#ifndef REPORT__H
#define REPORT__H 1

char *clocktime (char *atime, long seconds);
void writeline (Project * pr, char *s);
void writeheader (Project * pr, int type, int contin);
void writelogo (Project * pr);
int writereport (Project * pr);
void writesummary (Project * pr);
void writehydstat (Project * pr, int iter, double relerr);
void writemassbalance (Project * pr);
void writerelerr (Project * pr, int iter, double relerr);
void writestatchange (Project * pr, int k, char s1, char s2);
void writecontrolaction (Project * pr, int k, int i);
void writeruleaction (Project * pr, int k, char *ruleID);
int writehydwarn (Project * pr, int iter, double relerr);
void writehyderr (Project * pr, int errnode);
void writetime (Project * pr, char *fmt);

/**
 * function to clear contents of a project's report file.
 * \return error code.
 */
static inline int
clearreport (Project * pr)
{
  Report *rpt = &pr->report;
  if (rpt->RptFile == NULL)
    return 0;
  if (freopen (rpt->Rpt1Fname, "w", rpt->RptFile) == NULL)
    return 303;
  writelogo (pr);
  return 0;
}

/*
 * function to copy contents of a project's report file.
 *
 * \return error code.
 */
static inline int
copyreport (Report * rpt, 
            char *filename)     ///< name of file to copy to.
{
  FILE *tfile;
  int c;

  // Check that project's report file exists
  if (rpt->RptFile == NULL)
    return 0;

  // Open the new destination file
  tfile = fopen (filename, "w");
  if (tfile == NULL)
    return 303;

  // Re-open project's report file in read mode
  fclose (rpt->RptFile);
  rpt->RptFile = fopen (rpt->Rpt1Fname, "r");

  // Copy contents of project's report file
  if (rpt->RptFile)
    {
      while ((c = fgetc (rpt->RptFile)) != EOF)
        fputc (c, tfile);
      fclose (rpt->RptFile);
    }

  // Close destination file
  fclose (tfile);

  // Re-open project's report file in append mode
  rpt->RptFile = fopen (rpt->Rpt1Fname, "a");
  if (rpt->RptFile == NULL)
    return 303;
  return 0;
}

#endif
