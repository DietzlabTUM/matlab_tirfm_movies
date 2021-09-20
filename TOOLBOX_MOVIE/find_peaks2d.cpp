
 /*[filtered image] = my_highlowfilter[raw image, w, sigma, n ] n is number of standart defiation for gaussian filter  */
#include "mex.h"
#include "matrix.h"
#include "math.h"
#include "stdlib.h"
#include "peakprocessing.h"
#include <vector>


void mexFunction (
                  int nlhs,       mxArray *plhs[],
                  int nrhs, const mxArray *prhs[]
                  )
{
  
  int size_x, size_y, i,j,k,l, w,length, absolute;
  double *image,*max_image, *peak_out,  min_intensity, max_p;
  std::vector< std::vector<int> > peaks; 
  /* Check for proper number of arguments */
  if (nrhs != 4 ) {
    mexErrMsgTxt ("find_peaks requires an image, a width w (typically 3*sigma), a threshold, use absolute YES=1/NO=0");
  } else if (nlhs != 1 ) {
    mexErrMsgTxt ("add returns list of peaks ");
  }
 
  
  /* Get arguments from MATLAB */
  image = (mxGetPr(prhs[0]));
  size_x = mxGetM(prhs[0]);
  size_y = mxGetN(prhs[0]);
  w = (int)(mxGetScalar(prhs[1]));
  min_intensity = (mxGetScalar(prhs[2]));
  absolute = (int)(mxGetScalar(prhs[3]));

  PeakProcessing::find_peaks(image, 0, size_x, size_y, w, min_intensity, &peaks, absolute);

  std::vector< std::vector<double> > m0m2;
  PeakProcessing::noise_discrimination(image,0, size_x, size_y, 3*w, &peaks, &m0m2);

  length = (int)peaks.size();
	int rows = 0 ;
	if (peaks.size() > 0){
		rows = peaks[0].size();
	}
  /* Allocate storage place for the result */
  plhs[0] = mxCreateDoubleMatrix(length,rows, mxREAL); 
  peak_out = mxGetPr(plhs[0]);

  
//mexPrintf("find_peaks: dim = %i,%i \twidth = %i \th_min = %f", size_x, size_y, w, min_intensity);
if (absolute){
 //   mexPrintf("\tabsolute intensity\n", peaks.size());
           }
else {
//    mexPrintf("\trelative intensity\n", rows);
    }
// mexPrintf("Found %i peaks\n", peaks.size());

//fill the the outlist
  for(i=0; i< length; i++){
	  for (int j=0; j < rows; j++) {
		  peak_out[i+j*length] = peaks[i][j];
	  }
		/* peak_out[i+0*length] = peaks[i][0];
		 peak_out[i+1*length] = peaks[i][1]; 
		peak_out[i+2*length] = peaks[i][2];
		peak_out[i+3*length] = m0m2[i][0];
		peak_out[i+4*length] = m0m2[i][1]; */
		//mexPrintf("%i  %i  %i \n", peaks[i][0], peaks[i][1], peaks[i][2]);
  } 

//destroy peaks

  



}
