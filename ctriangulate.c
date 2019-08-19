#include <stdio.h>
#include <qhull_ra.h>
// This function is designed to be called from fortran 
// via the wrapper written in utils.F90  
//*-------------------------------------------------
//-call qhull to triangulate point set 
int ctriangulate(int DIM, int NUMPOINTS, double *fpoints) { // int *faces) {

   boolT ismalloc= False;    /* True if qhull should free points in qh_freeqhull() or reallocation */
   char flags[250];          /* option flags for qhull, see qh-quick.htm */
   FILE *outfile= NULL;      /* output from qh_produce_output()
                              use NULL to skip qh_produce_output() */
   FILE *errfile= stderr;    /* error messages from qhull code */
   int exitcode;             /* 0 if no error from qhull */
   int curlong, totlong;     /* memory remaining after qh_memfreeshort, used if !qh_NOmem  */
   int i;

   facetT *facet;            /* set by FORALLfacets */
   vertexT *vertex, **vertexp; 

    QHULL_LIB_CHECK
  {
    coordT points[DIM*NUMPOINTS]; /* array of coordinates for each point */

    // put points in the coordT type.
    // will have to check if this is absolutely necessary
    for ( i=0 ; i < NUMPOINTS; i++) {
       points[i*2]  =fpoints[i*2];
       points[i*2+1]=fpoints[i*2+1];
    }

    qhT qh_qh;                     /* Create an instance of Qhull */
    qhT *qh= &qh_qh;
    qh_zero(qh, errfile);

    sprintf(flags, "qhull d Tcv i");

    fflush(NULL);
    exitcode= qh_new_qhull(qh, DIM, NUMPOINTS, points, ismalloc,
                    flags, outfile, errfile);
    fflush(NULL);

    if (!exitcode)

      FORALLfacets { 
          if( !facet->upperdelaunay) {
            FOREACHvertex_(facet->vertices)
                printf(" %d", qh_pointid (qh, vertex->point) );
            printf("\n");
          }
      }

    qh_freeqhull(qh, !qh_ALL);                 /* free long memory */
    qh_memfreeshort(qh, &curlong, &totlong);  /* free short memory and memory allocator */
    if (curlong || totlong)
        fprintf(errfile, "qhull internal warnin: did not free %d bytes of long memory (%d pieces)\n", totlong, curlong);
    /* Exiting the block frees qh_qh */
  }

    
return exitcode;

}


