#include <stdio.h>
#include <stdlib.h>
#include <qhull_ra.h>

// These functions are called from FORTRAN via an interface in utils.F90


// This function is designed to be called from fortran 
// via the wrapper written in utils.F90  
//
//*-------------------------------------------------
//-call qhull to triangulate point set 
int *faces(int DIM, int NUMPOINTS, double *fpoints, int *NF) {

   boolT ismalloc= False;    /* True if qhull should free points in qh_freeqhull() or reallocation */
   char flags[250];          /* option flags for qhull, see qh-quick.htm */
   FILE *outfile= NULL;      /* output from qh_produce_output()
                              use NULL to skip qh_produce_output() */
   FILE *errfile= stderr;    /* error messages from qhull code */
   int exitcode;             /* 0 if no error from qhull */
   int curlong, totlong;     /* memory remaining after qh_memfreeshort, used if !qh_NOmem  */
   int nt=*NF; 
   int *tmp;
   int i,j;

   facetT *facet;            /* set by FORALLfacets */
   vertexT *vertex, **vertexp; 
  {
    qhT qh_qh;                     /* Create an instance of Qhull */
    qhT *qh= &qh_qh;
    qh_zero(qh, errfile);

    sprintf(flags, "qhull d Qt i");

    fflush(NULL);
    exitcode= qh_new_qhull(qh, DIM, NUMPOINTS, fpoints, ismalloc,
                    flags, outfile, errfile);
    fflush(NULL);

    if (!exitcode)
      tmp = (int*) malloc( sizeof(int)*(DIM+1)*(2*NUMPOINTS) );

      i=0;
      // we only need the lower hull facets (Delaunay ones)
      FORALLfacets { 
          if( !facet->upperdelaunay) {
            FOREACHvertex_(facet->vertices) {
                tmp[i]= qh_pointid (qh, vertex->point);
                i++;
                }
          }
      } 
      *NF = i/(DIM+1); 
      

    qh_freeqhull(qh, !qh_ALL);                 /* free long memory */
    qh_memfreeshort(qh, &curlong, &totlong);  /* free short memory and memory allocator */
    if (curlong || totlong)
        fprintf(errfile, "qhull internal warnin: did not free %d bytes of long memory (%d pieces)\n", totlong, curlong);
  }
// return facet table     
return tmp;

}

// 
void destroy_storage(int *ptr)
{
   free(ptr);
}


// This function is the pnpoly algorithm 
int pnpoly(int nvert, double *vertx, double *verty, double testx, double testy)
{
  int i, j, c = 0;
  for (i = 0, j = nvert-1; i < nvert; j = i++) {
    if ( ((verty[i]>testy) != (verty[j]>testy)) &&
	 (testx < (vertx[j]-vertx[i]) * (testy-verty[i]) / (verty[j]-verty[i]) + vertx[i]) )
       c = !c;
  }
  return c;
}


