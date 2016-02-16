/*
* Compute the discrete Frechet distance between two curves specified by
* discrete ordered points in n-dimensional space.
*
* Based on `DiscreteFrechetDist` by Zachary Danziger,
* http://www.mathworks.com/matlabcentral/fileexchange/ \
* 31922-discrete-frechet-distance
*
* This implementation omits the computation of the coupling sequence. Use
* this program if you only want to get the DFD, and get it fast.
*
* Implements algorithm from
* [1] T. Eiter and H. Mannila. Computing discrete Frechet distance.
*     Technical report 94/64, Christian Doppler Laboratory
*
* Copyright (c) 2016, Mikhail Pak
*/


/*
* Notice on compiler portability:
* Microsoft Windows SDK 7.1 seems to have problems with `fmin` and `fmax`.
* Compile this function using MinGW 4.9.2 C/C++ (TDM-GCC) by typing
* `mex DiscreteFrechetDistance.c`
*/


#include <math.h>
#include "mex.h"
#include "matrix.h"


/* Define functions */
/* Main function for computation of the discrete Frechet distance */
void discrete_frechet_distance(double *dfd);
/* Recursive function for computation of the coupling measure, see [1] */
double recursive_c(mwIndex i, mwIndex j);
/*
* Euclidean (2-) norm between the i-th point of the 1st curve and the j-th
* point of the 2nd curve; i = 1, ..., n_1; j = 1, ..., n_2.
*/
double norm_2(mwIndex i, mwIndex j);


/* Define global variables */
/*
* `c_1` and `c_2` : Arrays with the 1st and 2nd curve's points respectively
*
* Memory layout: Column-major (inherited from MATLAB), i.e. for `c_1`:
* [<p 1, d 1> .. <p 1, dim n_d> .. <p n_1, d 1> .. <p n_1, d n_d>]
* p -- point, d -- dimension
*/
double *c_1, *c_2;
/*
* `ca` : Search array (refer to [1], Table 1, matrix `ca`)
*
* Memory layout: Row-major
*/
double *ca;
/* `n_1` and `n_2` : Number of points of the 1st and 2nd curves */
mwIndex n_1, n_2;
/* `n_d` : Number of dimensions */
mwIndex n_d;


/* Entry point for the MEX interface */
void mexFunction(int nlhs, mxArray **plhs, int nrhs, const mxArray **prhs)
{
    double *dfd; /* Return value: discrete Frechet distance */
    const mwSize *size_c_1, *size_c_2; /* Size of the input arguments*/


    /* Target the pointer to the MEX output value (LHS) */
    plhs[0] = mxCreateDoubleScalar(mxREAL);
    dfd = mxGetPr(plhs[0]);
    /* Initialise it with a default value of -1.0 */
    *dfd = -1.0;

    /* Store the pointer to the input arrays with curve points*/
    c_1 = mxGetPr(prhs[0]);
    c_2 = mxGetPr(prhs[1]);

    /* Get size of the inputs*/
    size_c_1 = mxGetDimensions(prhs[0]);
    size_c_2 = mxGetDimensions(prhs[1]);

    /*
    * This function does not perform a lot of input checks. Please take
    * care of it in MATLAB.
    */

    /* Check only if dimensions of the two curves are consistent */
    if (size_c_1[0] != size_c_2[0]) {
        printf("Error: c_1 and c_2 must have an equal number of rows!\n");
        return;
    }

    /* Store array sizes to global variables */
    n_1 = size_c_1[1];
    n_2 = size_c_2[1];
    n_d = size_c_1[0];


    /* Call the main function */
    discrete_frechet_distance(dfd);
}


void discrete_frechet_distance(double *dfd)
{
    mwIndex k; /* Index for initialisation of `ca`*/

    /* Allocate memory for `ca` */
    ca = (double *) mxMalloc(n_1*n_2*sizeof(double));

    /* Initialise it with -1.0 */
    for (k = 0; k < n_1*n_2; k++)
    {
        *(ca + k) = -1.0;
    }

    /* Call the recursive computation of the coupling measure */
    *dfd = recursive_c(n_1, n_2);

    /* Free memory */
    mxFree(ca);
}


double recursive_c(mwIndex i, mwIndex j)
{
    double *ca_ij; /* Pointer to `ca(i, j)`, just to simplify notation */

    /*
    * Target the shortcut to the (i, j)-th entry of the matrix `ca`
    *
    * Once again, notice the 1-offset.
    */
    ca_ij = ca + (i - 1)*n_2 + (j - 1);

    /* This implements the algorithm from [1] */
    if (*ca_ij > -1.0)
    {
        return *ca_ij;
    }
    else if ((i == 1) && (j == 1))
    {
        *ca_ij = norm_2(1, 1);
    }
    else if ((i > 1) && (j == 1))
    {
        *ca_ij = fmax(recursive_c(i - 1, 1), norm_2(i, 1));
    }
    else if ((i == 1) && (j > 1))
    {
        *ca_ij = fmax(recursive_c(1, j - 1), norm_2(1, j));
    }
    else if ((i > 1) && (j > 1))
    {
        *ca_ij = fmax(
                      fmin(fmin(
                           recursive_c(i - 1, j    ),
                           recursive_c(i - 1, j - 1)),
                           recursive_c(i,     j - 1)),
                      norm_2(i, j));
    }
    else
    {
        /* Fetch the Inf value from the MEX library*/
        *ca_ij = mxGetInf();
    }

    return *ca_ij;
}


double norm_2(mwIndex i, mwIndex j)
{
    double dist, diff; /* Temp variables for simpler computations */
    mwIndex k; /* Index for iterating over dimensions */

    /* Initialise distance */
    dist = 0.0;

    for (k = 0; k < n_d; k++)
    {
        /*
        * Compute the distance between the k-th component of the i-th point
        * of the 1st curve and the k-th component of the j-th point of the
        * 2nd curve.
        *
        * Notice the 1-offset added for better readability (as in [1]).
        */
        diff = *(c_1 + (i - 1)*n_d + k) - *(c_2 + (j - 1)*n_d + k);
        /* Increment the accumulator variable with the squared distance */
        dist += diff*diff;
    }

    /* Compute the square root for the 2-norm */
    dist = sqrt(dist);

    return dist;
}
