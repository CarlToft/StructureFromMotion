
#include <math.h>
#include "mex.h"

/* Input Arguments */

#define	DESC1	prhs[0]
#define	DESC2	prhs[1]

/* Output Arguments */

#define	MATCHINGS	plhs[0]

static double inner(double *vec1, double *vec2, mwSize leng)
{
    int i;
    double sum;
    
    sum = 0;
    for (i=0; i < leng; i++)
    {
        sum += vec1[i]*vec2[i];
        /*sum = sum + vec1[i]*vec2[i];*/
        /*mexPrintf("%f %f\n",vec1[i],vec2[i]);*/
    }
    return sum;
}

static void mex_match(
		   double	*matchings,
		   double	*desc1,
           mwSize   rows1,
           mwSize   cols1,
 		   double	*desc2,
           mwSize   rows2,
           mwSize   cols2)
{
    int i,j,maxIndex,matchnr;
    double innerProd, maxInnerProd, sndInnerProd;
    
    matchnr = 0;
    for (i=0; i < cols1; i++) /* loopa igenom desc1 */
    {
        maxInnerProd = -1;
        sndInnerProd = -1;
        for (j = 0; j < cols2; j++)
        {
            innerProd = inner(desc1+i*rows1, desc2+j*rows2, rows1);
            
            /*mexPrintf("%f\n",innerProd);*/
            if (innerProd >= maxInnerProd)
            {
                sndInnerProd = maxInnerProd;
                maxInnerProd = innerProd;
                maxIndex = j+1;
            }
            else
            {
                if (innerProd > sndInnerProd)
                    sndInnerProd = innerProd;
            }
        }
        /*mexPrintf("%f %f\n",acos(sndInnerProd),acos(maxInnerProd));*/
        if (acos(maxInnerProd) < 0.5*acos(sndInnerProd))
        {
            matchings[i] = maxIndex;
            ++matchnr;
        }
    }
    mexPrintf("Found %d matches\n", matchnr);
    return;
}


void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray*prhs[] )
     
{ 
    double *desc1,*desc2,*matchings; 
    mwSize rows1,cols1,rows2,cols2; 
    
    /* Check for proper number of arguments */
    
    if (nrhs != 2) { 
	mexErrMsgTxt("Two input arguments required."); 
    } else if (nlhs > 1) {
	mexErrMsgTxt("Too many output arguments."); 
    } 
    
    rows1 = mxGetM(DESC1); 
    cols1 = mxGetN(DESC1);
    rows2 = mxGetM(DESC2); 
    cols2 = mxGetN(DESC2);
    
    /* Create a matrix for the return argument */ 
    MATCHINGS = mxCreateDoubleMatrix(1, cols1, mxREAL); 
    
    /* Assign pointers to the various parameters */ 
    matchings = mxGetPr(MATCHINGS);
    
    desc1 = mxGetPr(DESC1); 
    desc2 = mxGetPr(DESC2);
        
    /* Do the actual computations in a subroutine */
    
    mex_match(matchings,desc1,rows1,cols1,desc2,rows2,cols2); 
    return;
    
}


