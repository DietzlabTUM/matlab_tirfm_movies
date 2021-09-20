
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
	
	int col, row, r_find , movie_length, traces_col, min_length, trace_length, r_integrate, fnumber_size, trace_number;
	double *movie, h_min, *p2Matrix, *cur_trace, *fnumber, *pos;
	std::vector< std::vector< std::vector<int> > > itraces;
	std::vector< std::vector< int > > positions;
	
	std::vector<int> framenumber;
	/* Check for proper number of arguments */
	if (nrhs != 5 ) {
		mexErrMsgTxt ("append_traces_to_positions requires: movie, itraces, framenumbers, positions, r_integrate");
	} else if (nlhs != 1 ) {
		mexErrMsgTxt ("append_traces_to_positions returns list of itraces ");
	}
	
	
	/* Get arguments from MATLAB */
	//
	movie = mxGetPr(prhs[0]);  //columns * rows * frame
	movie_length = (int) mxGetDimensions(prhs[0])[2]; //length of video = number of frames
	col = (int) mxGetDimensions(prhs[0])[0]; //number of columns
	row = (int) mxGetDimensions(prhs[0])[1]; //number of rows
	
	r_integrate = (int)(mxGetScalar(prhs[4]));
	
	
	

	
	//fill framenumbers
	fnumber = mxGetPr(prhs[2]);  //columns * rows * frame
	for(int i=0; i< PeakProcessing::max(mxGetM(prhs[2]), mxGetN(prhs[2]) ); i++){
		framenumber.push_back( (int)fnumber[i] );
	}
	

	
	trace_number = mxGetM(prhs[1]);
	//fill itraces
	for(int i=0; i<trace_number; i++){
        cur_trace = mxGetPr(mxGetCell(prhs[1],i));
		trace_length = mxGetM(mxGetCell(prhs[1],i));
		std::vector< std::vector<int> > tmptrace;
		//loop through ith trace
		for(int j=0; j<trace_length; j++){
			std::vector<int> addthis;
			addthis.push_back(cur_trace[trace_length*0+j]); //framenumber
			addthis.push_back(cur_trace[trace_length*1+j]); //x
			addthis.push_back(cur_trace[trace_length*2+j]); //y
			addthis.push_back(cur_trace[trace_length*3+j]); //intensity
			tmptrace.push_back(addthis);
		} 
		itraces.push_back(tmptrace);
	}
	
	
	//load positions
	pos= mxGetPr(prhs[3]);  //columns * rows * frame
	int pos_col = mxGetM(prhs[3]);
	for(int i=0; i< pos_col; i++){
		std::vector< int > tmp;
		for (int j=0; j < mxGetN(prhs[3]) ; j++) {
			tmp.push_back( (int)pos[j*pos_col + i]);
		}
		
		positions.push_back(tmp);
	}
	
	
	PeakProcessing::append_traces_to_position(movie, movie_length, col, row,  &itraces, &positions, r_integrate, &framenumber);
	
	
	
	//delete short traces
	/*int trace_number=0;
	 for(int i=0; i< traces.size(); i++){
	 if(traces[i].size() < min_length){
	 traces.erase(traces.begin()+i);
	 i--;
	 }
	 } */
	
    mexPrintf("movie_length = %i \tImage = %i , %i\tr_integrate = %i\t%i itraces\n", movie_length, col, row, r_integrate, itraces.size());

	
	//Allocate storage place for the result 
	int dims[1];  dims[0]=(int) itraces.size();  //number of traces   
	plhs[0] = mxCreateCellArray(1, dims  );
	
	
	for(int i=0; i< itraces.size(); i++){
		traces_col = itraces[i].size();
		mxArray * tmp = mxCreateDoubleMatrix((int)itraces[i].size(),4, mxREAL); //make mxArray to store the ith trace
		p2Matrix = mxGetPr(tmp);
		//fill this array
		for(int j=0; j<itraces[i].size(); j++){
			p2Matrix[traces_col*0+j] = itraces[i][j][0]; //framenumber
			p2Matrix[traces_col*1+j] = itraces[i][j][1]; //x
			p2Matrix[traces_col*2+j] = itraces[i][j][2]; //y 
			p2Matrix[traces_col*3+j] = itraces[i][j][3]; //intensity 

			//mexPrintf("%i\t %i \t %i\n", traces[i][j][0], traces[i][j][1], traces[i][j][2]);
		}
		mxSetCell(plhs[0], i,tmp);
		
	}
	
	
	
	
	
	
	
}
