//
//  Created by Jean-Marie Mirebeau on 26/06/2015.
//  Copyright (c) 2015 Jean-Marie Mirebeau. All rights reserved.
//

//#include <cassert>
#include "mex.h"

size_t Index(size_t i,size_t j,size_t rows){return i+j*rows;}

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[] ){
    if(nrhs!=3)
        mexErrMsgIdAndTxt( "TripletsProduct:invalidNumInputs",
                          "Three arguments needed. (triplets, vector, [target time,safety]).");

    for(int i=0; i<3; ++i)
        if(!mxIsNumeric(prhs[i]) || !mxIsDouble(prhs[i]))
            mexErrMsgIdAndTxt( "TripletsProduct:invalidnputs",
                              "Inputs must be numeric arrays of double type.");
    
    const mxArray * mxTriplets = prhs[0], * mxVector = prhs[1], * mxTimeData = prhs[2];
    const double * pTriplets = mxGetPr(mxTriplets), * pVector=mxGetPr(mxVector), * pTimeData = mxGetPr(mxTimeData);
    
    const mwSize rows = mxGetM(mxVector), cols = mxGetN(mxVector), nTriplets = mxGetN(mxTriplets);
    if(mxGetM(mxTriplets)!=3)
        mexErrMsgIdAndTxt("TripletsProduct:invalidInputs","Invalid triplets size");
    
    if(nlhs!=2)
        mexErrMsgIdAndTxt("TripletsProduct:invalidOutputs",
                          "Two ouputs needed. (vector, time reached)");
    
    plhs[0] = mxCreateDoubleMatrix( rows, cols, mxREAL);
    plhs[1] = mxCreateDoubleMatrix( 1, 1, mxREAL);
    
    double *pResult = mxGetPr(plhs[0]), *reached = mxGetPr(plhs[1]);;
    for(size_t k=0; k<rows*cols; ++k) pResult[k]=0.; *reached=0; //std::fill(pResult, pResult+rows*cols, 0.);
    
    for(size_t k=0; k<nTriplets; ++k){
        const double i_ = -1+pTriplets[Index(0,k,3)], j_ = -1+pTriplets[Index(1,k,3)], c = pTriplets[Index(2,k,3)];
        const size_t i=size_t(i_), j=size_t(j_);
        
        if(double(i)!=i_ || double(j)!=j_)      mexErrMsgIdAndTxt("TripletsProduct:tripletsFormat","TripletsProduct:indices are not integers");
        if(!(i<rows && j<rows)) mexErrMsgIdAndTxt("TripletsProduct:tripletsFormat","TripletsProduct:indices are out of range");
        
        for(int k=0; k<cols; ++k)
            pResult[Index(i,k,cols)] += c*pVector[Index(j,k,cols)];
    }
    
}