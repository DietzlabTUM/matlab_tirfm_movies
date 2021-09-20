
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
	double *movie, h_min, *p2Matrix, *cur_trace, *fnumber;
	std::vector< std::vector< std::vector<int> > > traces;
	std::vector< std::vector< std::vector<int> > > itraces;

	std::vector<int> framenumber;
	/* Check for proper number of arguments */
	if (nrhs != 8 ) {
		mexErrMsgTxt ("append_traces requires: movie, traces, itraces, framenumbers, h_min, r_find, r_integrate, min_length");
	} else if (nlhs != 2 ) {
		mexErrMsgTxt ("add returns list of traces ");
	}
	
	
	/* Get arguments from MATLAB */
	//
	movie = mxGetPr(prhs[0]);  //columns * rows * frame
	movie_length = (int) mxGetDimensions(prhs[0])[2]; //length of video = number of frames
	col = (int) mxGetDimensions(prhs[0])[0]; //number of columns
	row = (int) mxGetDimensions(prhs[0])[1]; //number of rows
	h_min = (mxGetScalar(prhs[4]));
	r_find = (int)(mxGetScalar(prhs[5]));
	r_integrate = (int)(mxGetScalar(prhs[6]));
	min_length = (int)(mxGetScalar(prhs[7]));


	
		//fill framenumbers
	fnumber = mxGetPr(prhs[3]);  //columns * rows * frame
	for(int i=0; i< PeakProcessing::max(mxGetM(prhs[3]), mxGetN(prhs[3]) ); i++){
		framenumber.push_back((int)fnumber[i]);
	}
	
	trace_number = mxGetM(prhs[1]);
	//fill traces
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
			tmptrace.push_back(addthis);
		} 
		traces.push_back(tmptrace);
	}
	
	trace_number = mxGetM(prhs[2]);
	//fill itraces
	for(int i=0; i<trace_number; i++){
        cur_trace = mxGetPr(mxGetCell(prhs[2],i));
		trace_length = mxGetM(mxGetCell(prhs[2],i));
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
	
	
	

	PeakProcessing::append_traces(movie, movie_length, col, row, &traces, &itraces, r_find, h_min, r_integrate, &framenumber);
	

	
	
	//delete short traces
	for(int i=0; i< traces.size(); i++){
		if(traces[i].size() < min_length){
			traces.erase(traces.begin()+i);
			itraces.erase(itraces.begin()+i);
			i--;
		}
	} 
	
    mexPrintf("movie_length = %i \tImage = %i , %i \t", movie_length, col, row);
	mexPrintf("r_find = %i \t h_min = %.2f \t r_integrate = %i\t", r_find, h_min, r_integrate);
	mexPrintf("%i traces\n", traces.size());

	//Allocate storage place for the result 
	int dims[1];  dims[0]=(int) traces.size();  //number of traces   
	plhs[0] = mxCreateCellArray(1, dims  );
	
	
	for(int i=0; i< traces.size(); i++){
		traces_col = traces[i].size();
		mxArray * tmp = mxCreateDoubleMatrix((int)traces[i].size(),3, mxREAL); //make mxArray to store the ith trace
		p2Matrix = mxGetPr(tmp);
		//fill this array
		for(int j=0; j<traces[i].size(); j++){
			p2Matrix[traces_col*0+j] = traces[i][j][0]; //framenumber
			p2Matrix[traces_col*1+j] = traces[i][j][1]; //x
			p2Matrix[traces_col*2+j] = traces[i][j][2]; //y 
			//mexPrintf("%i\t %i \t %i\n", traces[i][j][0], traces[i][j][1], traces[i][j][2]);
		}
		mxSetCell(plhs[0], i,tmp);
		
	}
	
	
	
	//Allocate storage place for the result 
	dims[0]=(int) itraces.size();  //number of traces   
	plhs[1] = mxCreateCellArray(1, dims  );
	
	
	for(int i=0; i< itraces.size(); i++){
		traces_col = itraces[i].size();
		mxArray * tmp = mxCreateDoubleMatrix((int)itraces[i].size(),4, mxREAL); //make mxArray to store the ith trace
		p2Matrix = mxGetPr(tmp);
		//fill this array
		for(int j=0; j<itraces[i].size(); j++){
			p2Matrix[traces_col*0+j] = itraces[i][j][0]; //framenumber
			p2Matrix[traces_col*1+j] = itraces[i][j][1]; //x
			p2Matrix[traces_col*2+j] = itraces[i][j][2]; //y 
			p2Matrix[traces_col*3+j] = itraces[i][j][3]; //y 

			//mexPrintf("%i\t %i \t %i\n", traces[i][j][0], traces[i][j][1], traces[i][j][2]);
		}
		mxSetCell(plhs[1], i,tmp);
		
	}
	
	
	
	
	
}
