

#ifndef PEAKPROCESSING_H
#define PEAKPROCESSING_H
#include <vector>

class PeakProcessing
{
    public:
	//finding peaks
	static void find_peaks( double * , int,int,int, int, double, std::vector<std::vector<int> >  *, int absolute); // video,offset, size_x, size_y, width, threshold, peaks_out, use absolute height
	static void refine_peaks(double * , int,int , int , int ,std::vector<std::vector<int> >  * );  //peaks, width w   <--------- if one needs accurate central position... i  prob wont need this
	static void noise_discrimination(double *, int,int, int, int, std::vector<std::vector<int> >  *, std::vector<std::vector<double> >  *); //image, size_x, size_y, width w, (filled)peaks, m2-m0-out

	//connecting peaks
	static void makeTraces(double *, int, int, int, int, int, int , double, std::vector< std::vector< std::vector<int> > > *); //video, size_x, size_y, length, start, skip, width, treshold, tracelist_out
	static void mergeTraces(std::vector< std::vector< std::vector<int> > > *, std::vector< std::vector< std::vector<int> > > *,std::vector< std::vector< std::vector<int> > > *,double , std::vector< std::vector< std::vector< std::vector<int> > > > *);//trace1, trace2, ,mergedtrace

	static void merge2Traces(std::vector< std::vector< std::vector<int> > > *, std::vector< std::vector< std::vector<int> > > *, double , std::vector< std::vector< std::vector< std::vector<int> > > > *);//trace1, trace2, ,mergedtrace

	static void integrate(double *, int, int, int, int, int, int,  std::vector< std::vector<int> >  *,  std::vector< std::vector<int> >  * ); //movie, movie_length, size_x, size_y, width, start, skip, trace, trace_out
    
	static void integrate_double(double *, int, int, int, int, int, int,  std::vector< std::vector<double> >  *,  std::vector< std::vector<double> >  * ); //movie, movie_length, size_x, size_y, width, start, skip, trace, trace_out
   
	static void append_traces(double * movie,int movie_length,int size_x, int size_y, std::vector< std::vector< std::vector<int> > > * traces, std::vector< std::vector< std::vector<int> > > * itraces, double r_find, double h_min, double r_integrate, std::vector<int> * framenumber);

	static void append_traces_to_position(double * movie,int movie_length,int size_x, int size_y, std::vector< std::vector< std::vector<int> > > * itraces, std::vector< std::vector<int> > * positions, double r_integrate, std::vector<int> * framenumber);

	static void map_traces(std::vector< std::vector<double> > * pos1, std::vector< std::vector<double> > * pos2, std::vector< std::vector<double> > * pos3, std::vector< std::vector<int> > * map, double r);

	
	static int max(int ,int);
	static int min(int ,int);
	
	private:
	static void connect(std::vector< std::vector< std::vector<int> > > * , std::vector< std::vector<int> > * , double, int); //traces, peaks in image: links the peaks in peaks to the traces in traces
	static double mylength(int, int);
	
	static void average_position(std::vector< std::vector< std::vector<int> > > * , std::vector< std::vector<double> >* );
	
	
};

#endif
