
/*[traces] = make_traces[video, start, skip, width,  ] n is number of standart defiation for gaussian filter  */
#include "mex.h"
#include "matrix.h"
#include "math.h"
#include "stdlib.h"
#include "peakprocessing.h"
#include <vector>


void mexFunction (
                  int nlhs,       mxArray *plhs[],
                  int nrhs,  const mxArray *prhs[]
                  )
{
	
	int pos_col;
	double  *pos1, *pos2,*pos3, *map_out, r;
	std::vector< std::vector< double > > positions1;
	std::vector< std::vector< double > > positions2;
	std::vector< std::vector< double > > positions3;
	std::vector< std::vector< int > > map;

	
	std::vector<int> framenumber;
	/* Check for proper number of arguments */
	if (nrhs != 4 ) {
		mexErrMsgTxt ("map_traces requires: pos1 (N*2),pos2 (N*2), pos3 (N*2), radius ");
	} else if (nlhs != 1 ) {
		mexErrMsgTxt ("map_traces returns a mapping matrix ");
	}
	
	r = (int)(mxGetScalar(prhs[3]));

	
	
	//load positions
	pos1= mxGetPr(prhs[0]);  //columns * rows * frame
	pos_col = mxGetM(prhs[0]);
	for(int i=0; i< pos_col; i++){
		std::vector< double > tmp;
		for (int j=0; j < mxGetN(prhs[0]) ; j++) {
			tmp.push_back( pos1[j*pos_col + i]);
		}
		positions1.push_back(tmp);
	}
	
	//load positions
	pos2= mxGetPr(prhs[1]);  //columns * rows * frame
	pos_col = mxGetM(prhs[1]);
	for(int i=0; i< pos_col; i++){
		std::vector< double > tmp;
		for (int j=0; j < mxGetN(prhs[1]) ; j++) {
			tmp.push_back(pos2[j*pos_col + i]);
		}
		positions2.push_back(tmp);
	}
	//load positions
	pos3= mxGetPr(prhs[2]);  //columns * rows * frame
	pos_col = mxGetM(prhs[2]);
	for(int i=0; i< pos_col; i++){
		std::vector< double > tmp;
		for (int j=0; j < mxGetN(prhs[2]) ; j++) {
			tmp.push_back(pos3[j*pos_col + i]);
		}
		positions3.push_back(tmp);
	}
	
	/*
	mexPrintf("r = %.2f\n", r);
	
	
	mexPrintf("Positions 1:\n");
	for (int i=0; i<positions1.size(); i++) {
		mexPrintf("%.2f\t%.2f\n", positions1[i][0], positions1[i][1]);
	}
	mexPrintf("Positions 2:\n");

	for (int i=0; i<positions2.size(); i++) {
		mexPrintf("%.2f\t%.2f\n", positions2[i][0], positions2[i][1]);
	}
	mexPrintf("Positions 3:\n");

	for (int i=0; i<positions3.size(); i++) {
		mexPrintf("%.2f\t%.2f\n", positions3[i][0], positions3[i][1]);
	}
	 */
	
	
	
	
	
	PeakProcessing::map_traces(&positions1, &positions2, &positions3, &map, r);
	
	
		mexPrintf("%i,%i,%i merged to %i\n", positions1.size(),positions2.size(),positions3.size(), map.size());

	

	
	//make the output	 
	plhs[0] = mxCreateDoubleMatrix(map.size(),3, mxREAL); 
  	map_out = mxGetPr(plhs[0]);
	
	
	//fill output
	for(int i=0; i<map.size(); i++){
		for (int j=0; j < map[i].size(); j++) {
			map_out[j*map.size()+i] = map[i][j];
		}
		
	}
	
	
	
	
	
	
	
}
