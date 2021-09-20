#include "peakprocessing.h"
#include "math.h"
#include "stdlib.h"
#include "mex.h"
#include <vector>
#include <set>

using namespace std;

int PeakProcessing::max(int a, int b){
	if(b>a) return b;
	else return a;
}

int PeakProcessing::min(int a, int b){
	if(b<a) return b;
	else return a;
}

double PeakProcessing::mylength(int a, int b){
	return sqrt(pow(a,2) + pow(b,2));
}



void PeakProcessing::find_peaks( double * video,int offset, int size_x, int size_y, int w, double min_intensity,  std::vector< std::vector<int> > * peaks, int absolute){
	int i,j, w_bg, counter, wk, k_found, l_found;
	double local_max, global_min, bg;
	w_bg = 5*w;

	double n;
	double lin, var, I, I_bg;
	
    global_min = video[offset]; // init threshold to first pixel
    
	bool is_candidate[size_x][size_y]; // this might save some time as soon as a local max was found the algorithm will hopp to the next possible location
	for(i=0; i<size_x; i++){
		for(j=0; j<size_y; j++){
			is_candidate[i][j] = true;
            if (video[i*size_y+j+offset] < global_min) { // find global minimum
                global_min = video[i*size_y+j+offset];
            }
		}
	}

	for(i=0; i < size_x; i++){
	for(j=0; j < size_y; j++){
		//for each pixel (i,j)
		local_max = global_min; // local maximal value
		
		if (is_candidate[i][j]) { // if the spot is a potential candidate
			//loop through vicinity of pixel(i,j) and look for the maximum
			for(int k=max(0,i-w); k<=min(size_x-1,i+w);k++){ 	
				for(int l=max(0,j-w); l<=min(size_y-1,j+w);l++){
					if( video[k*size_y+l+offset] > local_max){
						local_max = video[k*size_y+l+offset];
                        k_found = k; l_found = l;
					}
				}
			} 
			
			
			
			
			if( (local_max == video[i*size_y+j+offset]) && (k_found == i) && (l_found == j) ){  //add pixel to peaklist if it is a local maxima and larger than min threshold
				
				//remove candidates, points in the vicinity are not local maxima
				for(int k=max(0,i-w); k<=min(size_x-1,i+w);k++){ 	
					for(int l=max(0,j-w); l<=min(size_y-1,j+w);l++){
						is_candidate[k][l] = false; //points in vicinity are not candidates
					}
				} 
				
				
				
				//subtract background from local maxima i,j&& (local_max >= min_intensity), calculate background intensity
				bg = 0.;
				I_bg = 0; counter = 0;
				for(int k=max(0,i-w_bg); k<=min(size_x-1,i+w_bg);k++){ 	
					for(int l=max(0,j-w_bg); l<=min(size_y-1,j+w_bg);l++){
						bg += video[k*size_y+l+offset];	
						I_bg += video[k*size_y+l+offset];
						counter++;
					}
				} 				
				bg /= counter; // /= ((double)(min(size_x-1,i+w_bg)-max(0,i-w_bg)+1)*(min(size_y-1,j+w_bg)-max(0,j-w_bg)+1));//normalize bg
				I_bg /= counter;
				//calculate intensity of the spot
				I = 0; counter = 0;
				for(int k=max(0,i-w); k<=min(size_x-1,i+w);k++){ 	
					wk = (int) (  sqrt(pow(w,2)-pow(i-k,2))  ); //try circular
					for(int l=max(0,j-wk); l<=min(size_y-1,j+wk);l++){
						I += (int) (video[k*size_y+l+offset]);
						counter++;
					}
				} 
				I_bg = I_bg * counter;
				
				/*
				n = ((double)(min(size_x-1,i+w_bg)-max(0,i-w_bg)+1)*(min(size_y-1,j+w_bg)-max(0,j-w_bg)+1));
				for(int k=max(0,i-w_bg); k<=min(size_x-1,i+w_bg);k++){ 	
					for(int l=max(0,j-w_bg); l<=min(size_y-1,j+w_bg);l++){
						lin = fabs(bg - video[k*size_y+l+offset]);
						var = pow(bg - video[k*size_y+l+offset], 2);
					}
				} 
				lin = lin/n;
				var = var/(n-1); */
				
				
				// if abolute == true => compare absolute intensity, else => compare relative intesity
				if((local_max-bg >= min_intensity && !absolute) || (local_max >= min_intensity && absolute) ){
					//mexPrintf("added: %i %i %.2f\n", i, j, local_max-bg );
					std::vector<int> tmp; 
					tmp.push_back(i);
					tmp.push_back(j);
					tmp.push_back(video[i*size_y+j+offset]);
					tmp.push_back(video[i*size_y+j+offset]-bg);
					tmp.push_back(I);
					tmp.push_back(I-I_bg);
					peaks->push_back(tmp); 
					//mexPrintf("%i,%i\tbg=%.2f\tlin=%.2f\tvar=%.2f\tsigma=%.2f\tmin_h=%.2f\tI=%.2f\n",i,j,bg, lin, var, sqrt(var), min_intensity, local_max-bg);
				}
			}
		}
		 
		
	}
  	}

	//refine_peaks(video, offset, size_x, size_y, w, peaks);
}


void PeakProcessing::makeTraces(double * video, int col, int row ,int length, int start, int skip, int w, double min_intensity, std::vector< std::vector< std::vector<int> > > * traces){ //video, size_x, size_y, start, skip, width, treshold, tracelist
	int i;
	vector< vector<int> > peaks;
	vector< vector<int> > tmp;
	if(start < length){
		i=start;
	}
	else{
		mexPrintf("Length of video is smaller than first index. start=%i\n", 0);
		i=0;
	}
	//mexPrintf("Procssing frame... ");
	while(i < length){ //i index of video
		tmp.clear(); //clear the last peaks vector
		//mexPrintf("%i ", i);
 		
							

		PeakProcessing::find_peaks(video, col*row*i, col, row , w, min_intensity, &tmp, 0);	
		//print peaks
		/* mexPrintf("Frame %i:\t", i);
		for(int j=0; j<tmp.size(); j++){
			mexPrintf("(%i %i) , ",tmp[j][0],tmp[j][1]);
		}  */ 

		

		PeakProcessing::connect(traces, &tmp, (double) w, i); //check if w works here
		
		/*mexPrintf("\nTraces:\n");
		for(int j=0; j<traces->size(); j++){
			mexPrintf("Trace: %i\n", j);
 			for(int k=0; k<traces->at(j).size(); k++){
				mexPrintf("%i\t%i\t%i\n",traces->at(j).at(k).at(0), traces->at(j).at(k).at(1), traces->at(j).at(k).at(2));
			}
		} */


		
		i += skip+1; //skip some frames
	}

}


void PeakProcessing::connect(std::vector< std::vector< std::vector<int> > > * traces, std::vector< std::vector<int> > * peaks, double d, int framenumber){	
	int j,x,y,dx,dy,last, j_min;	
	double d_min;
	//mexPrintf("d is %.2f\n", d);

	vector<int> free_traces; //containes free indexes of traces
	for(int i=0; i<traces->size(); i++){
		free_traces.push_back(i);
	}

	//loop through current peaks
	while(peaks->size() > 0){
		x = peaks->back().at(0);  //look at last element in current peak list
		y = peaks->back().at(1);
		d_min = d+1.;
		j_min = -1;
				
		//look for the closest previous peak in traces ... loop over free traces
		for(int i=0; i < (int) free_traces.size(); i++){ 
			j = free_traces[i]; //trace index
			/*
			last = traces->at(j).size()-1;
			dx = x - traces->at(j).at(last).at(1);
			dy = y - traces->at(j).at(last).at(2);
			if( ( mylength(dx,dy) <= d ) && ( mylength(dx,dy) < d_min )  ){
				j_min = j;
				d_min = mylength(dx,dy) ;
			} */
			for(int k=0; k<traces->at(j).size(); k++){
				dx = x - traces->at(j).at(k).at(1);
				dy = y - traces->at(j).at(k).at(2);
				if( ( mylength(dx,dy) <= d ) && ( mylength(dx,dy) < d_min )  ){
					j_min = j;
					d_min = mylength(dx,dy) ;
				}

			}
			
		}
				
		//current peak to be added
		vector<int> addthis;
		addthis.push_back(framenumber); //framenumber
		addthis.push_back(x); //x coordinate of peak 
		addthis.push_back(y); //y coordinate of peak 
		//mexPrintf("%i\t%i\t%i\n", framenumber, x, y);

		if(j_min < 0 ){ //there is no matching previous peak... start a new trace
			vector< vector<int> > newtrace;
			newtrace.push_back(addthis);
			traces->push_back(newtrace);
			j_min =  traces->size()-1; //remember in order to remove from free_traces		
		}
		else{ //append current peak to the trace 
			traces->at(j_min).push_back(addthis);  //add to trace i
		}
		
		//remove j_min from free traces
		for(int i=0; i < (int) free_traces.size(); i++){ 
			if(free_traces[i] == j_min){
				free_traces.erase(free_traces.begin()+i);
			}	
			
		}
		peaks->pop_back(); //remove last element in current peak list

		//problem?: two peaks in prev are close and one dissapears... there might be a wrong aiggnement <------------------- might be a PROBLEM
	}	

	
}

void PeakProcessing::noise_discrimination(double * image, int offset, int size_x, int size_y, int w, std::vector<std::vector<int> >  * peaks, std::vector<std::vector<double> >  * m2m0){ //image, size_x, size_y, width w, (filled)peaks, m2-m0-out
	int i,k,l, x,y, wk, counter;
	double m0, m2;
	for(i=0; i<peaks->size(); i++){ //loop throug peaks
		x = peaks->at(i).at(0);
		y = peaks->at(i).at(1);
		m0=0;
		m2=0;
		counter=0;
		for(k=max(0,x-w); k<=min(size_x-1,x+w);k++){ 	
			wk = (int) (  sqrt(pow(w,2)-pow(x-k,2))); //try circular
			for(l=max(0,y-wk); l<=min(size_y-1,y+wk);l++){
				m0 += image[offset+l+k*size_y];
				m2 += (image[offset+l+k*size_y])*(pow(x-k,2)+pow(y-l,2));
				counter++;
			}
		} 
		
		m0 = m0 / ((double) counter);
		if(fabs(m0)<1){
			m2 = 0;
		}
		else{
			m2 = m2/m0 ;
		}
		/*if(fabs(m2)>500){
			//nothing here
			m2=0.;
		} */
		
		vector<double> tmp;
		tmp.push_back(m0);
		tmp.push_back(m2);		
		m2m0->push_back(tmp);
	}
}



void PeakProcessing::refine_peaks(double * image, int offset, int size_x, int size_y, int w,std::vector<std::vector<int> >  * peaks){  //peaks, width w   <--------- if one needs accurate central position... i  prob wont need this
	int i,k,l, x,y, wk;
	double m0, mx, my;
	for(i=0; i<peaks->size(); i++){ //loop throug peaks
		x = peaks->at(i).at(0);
		y = peaks->at(i).at(1);
		m0=0; 
		mx=0;
		my=0;
		for(k=max(0,x-w); k<=min(size_x-1,x+w);k++){ 	
			
			wk = (int) (  sqrt(pow(w,2)-pow(x-k,2))  ); //try circular
			
			for(l=max(0,y-wk); l<=min(size_y-1,y+wk);l++){
				m0 += (double) image[l+k*size_y+offset];
				mx += image[l+k*size_y+offset] * (x-k);
				my += image[l+k*size_y+offset] * (y-l);
			}
			
		} 

		
		if(m0 > 0){
			peaks->at(i).at(0) += (int) round(mx / m0);
			peaks->at(i).at(1) += (int) round(my / m0);
		}
	
		/* if((int) round(mx / m0) > 0   || (int) round(my / m0) > 0 ){
			mexPrintf("%i   %i   dx = %i       dy = %i\n",x, y, (int) round(mx / m0), (int) round(my / m0) );
		} */
	}

}


/* this function takes 3 vectors of traces: trace1, trace2, trace3 and merges them into a merged_traces vector */
void PeakProcessing::mergeTraces(std::vector< std::vector< std::vector<int> > > * trace1, std::vector< std::vector< std::vector<int> > > *trace2, std::vector< std::vector< std::vector<int> > > *trace3, double w, std::vector< std::vector< std::vector< std::vector<int> > > > *merged_traces){//trace1, trace2, ,mergedtrace
	

	double d12,d13,d23, dx12, dy12, dx13, dy13, dx23, dy23, s1x2, s1y2, s2x2, s2y2, s3x2, s3y2, d_min;
	int j_min, k_min;
	vector< vector<double> > avg_1;
	vector< vector<double> > avg_2;
	vector< vector<double> > avg_3;
	PeakProcessing::average_position(trace1, &avg_1);
	PeakProcessing::average_position(trace2, &avg_2);
	PeakProcessing::average_position(trace3, &avg_3);
	merged_traces->clear();
	
	
	//find corresponding traces
	j_min = -1;
	k_min = -1;
	for(int i=0; i<avg_1.size(); i++){ //loop through avg1
		s1x2=pow(avg_1[i][8],2); 
		s1y2=pow(avg_1[i][9],2);
		j_min =-1;k_min =-1;
		d_min = sqrt(3)*w +1.;
		for(int j=0; j<avg_2.size(); j++){ //loop through avg2
			dx12=pow(avg_1[i][3]-avg_2[j][1], 2);
			dy12=pow(avg_1[i][4]-avg_2[j][2], 2);
			d12 = dx12+dy12;
			s2x2 = pow(avg_2[j][4],2);
			s2y2 = pow(avg_2[j][5],2);
			
			if( sqrt(d12) < w){ //point 1 and point 2 are at least w apart


			for(int k=0; k<avg_3.size(); k++){
				s3x2 = pow(avg_3[k][4],2);
				s3y2 = pow(avg_3[k][5],2);
				dx13 = pow(avg_1[i][3]-avg_3[k][1], 2);
				dy13 = pow(avg_1[i][4]-avg_3[k][2], 2);
				dx23 = pow(avg_2[j][1]-avg_3[k][1], 2);
				dy23 = pow(avg_2[j][2]-avg_3[k][2], 2);
				d23 = dx23+dy23;
				d13 = dx13+dy13;
				if(sqrt(d23)<w && sqrt(d13)<w){
					if(sqrt(d12+d13+d23) < d_min){
						j_min = j;
						k_min = k;
						d_min = sqrt(d12+d13+d23);
					}
				}
			}
			
	

			}
		}

		if(j_min >= 0 && k_min >=0){ //combine trace_1[i] trace_2[j_min] trace3[k_min]
			vector< vector< vector<int> > >	tmp_traces;
			vector< vector<int> > green;
			vector< vector<int> > red;
			vector< vector<int> > fret;
			//green		
			for(int j=0; j<trace1->at(i).size(); j++){
				vector<int> addthis;				
				for(int k=0; k< trace1->at(i).at(j).size(); k++){
					addthis.push_back(trace1->at(i).at(j).at(k));
				}
				green.push_back(addthis);
			}

			//red		
			for(int j=0; j<trace2->at(j_min).size(); j++){
				vector<int> addthis;				
				for(int k=0; k< trace2->at(j_min).at(j).size(); k++){
					addthis.push_back(trace2->at(j_min).at(j).at(k));
				}
				red.push_back(addthis);
			}

			//fret		
			for(int j=0; j<trace3->at(k_min).size(); j++){
				vector<int> addthis;
				for(int k=0; k< trace3->at(k_min).at(j).size(); k++){
					addthis.push_back(trace3->at(k_min).at(j).at(k));
				}
				fret.push_back(addthis);
			}
			tmp_traces.push_back(green);
			tmp_traces.push_back(red);
			tmp_traces.push_back(fret);
			merged_traces->push_back(tmp_traces);
		}
		//else there is no trace to add

	}



}


void PeakProcessing::merge2Traces(std::vector< std::vector< std::vector<int> > > * trace1, std::vector< std::vector< std::vector<int> > > *trace2, double w, std::vector< std::vector< std::vector< std::vector<int> > > > *merged_traces){
	
	double d12,d13,d23, dx12, dy12, dx13, dy13, dx23, dy23, s1x2, s1y2, s2x2, s2y2, s3x2, s3y2, d_min;
	int j_min, k_min;
	vector< vector<double> > avg_1;
	vector< vector<double> > avg_2;
	PeakProcessing::average_position(trace1, &avg_1);
	PeakProcessing::average_position(trace2, &avg_2);
	merged_traces->clear();
	
	
	//find corresponding traces
	j_min = -1;
	k_min = -1;
	for(int i=0; i<avg_1.size(); i++){ //loop through avg1
		s1x2=pow(avg_1[i][8],2); 
		s1y2=pow(avg_1[i][9],2);
		j_min =-1;
		d_min = sqrt(3)*w +1.;
		for(int j=0; j<avg_2.size(); j++){ //loop through avg2
			dx12=pow(avg_1[i][3]-avg_2[j][1], 2);
			dy12=pow(avg_1[i][4]-avg_2[j][2], 2);
			d12 = dx12+dy12;
			s2x2 = pow(avg_2[j][4],2);
			s2y2 = pow(avg_2[j][5],2);
			
			if( sqrt(d12) < w){ //point 1 and point 2 are at least w apart
				j_min = j;
				d_min = sqrt(d12);
			}
		}

		if(j_min >= 0 ){ //combine trace_1[i]  and trace_2[j_min] 
			vector< vector< vector<int> > >	tmp_traces;
			vector< vector<int> > green;
			vector< vector<int> > red;
			//green		
			for(int j=0; j<trace1->at(i).size(); j++){
				vector<int> addthis;				
				for(int k=0; k< trace1->at(i).at(j).size(); k++){
					addthis.push_back(trace1->at(i).at(j).at(k));
				}
				green.push_back(addthis);
			}

			//red		
			for(int j=0; j<trace2->at(j_min).size(); j++){
				vector<int> addthis;				
				for(int k=0; k< trace2->at(j_min).at(j).size(); k++){
					addthis.push_back(trace2->at(j_min).at(j).at(k));
				}
				red.push_back(addthis);
			}

			
			tmp_traces.push_back(green);
			tmp_traces.push_back(red);
			merged_traces->push_back(tmp_traces);
		}
		//else there is no trace to add

	}



}


void PeakProcessing::average_position(std::vector< std::vector< std::vector<int> > > * traces, std::vector< std::vector<double> > * out){
	vector<double> x;
	vector<double> x2;
	double n;
	out->clear();
	for(int i=0; i < traces->size() ; i++){
		x.clear();
		x2.clear();
		vector<double> addthis;
		n = (double) traces->at(i).size();
		

			for(int j=0; j<traces->at(i).at(0).size(); j++){
			x.push_back(0.);
			x2.push_back(0.);
			}
			
			for(int j=0; j< traces->at(i).size(); j++){
			for(int k=0; k<traces->at(i).at(0).size(); k++){
				x[k] += traces->at(i).at(j).at(k);
				x2[k] += pow(traces->at(i).at(j).at(k), 2);
			}			
			}
		
		
		
		for(int j=0; j<traces->at(i).at(0).size(); j++){
			addthis.push_back(x[j]/n); //append average values
		}
		if (n > 1) {
			for(int j=0; j<traces->at(i).at(0).size(); j++){
				addthis.push_back(  sqrt( x2[j]-pow(x[j],2)/n ) / (n-1) );  //append sigma
			}
		}
		else {
			for(int j=0; j<traces->at(i).at(0).size(); j++){
				addthis.push_back(  0.0 );  //append sigma
			}
		}

		

		
		
		out->push_back(addthis);
	}
}



void PeakProcessing::integrate(double * movie, int movie_length, int size_x, int size_y, int w, int start, int skip,  std::vector< std::vector<int> >  * trace,  std::vector< std::vector<int> >  * intensity_trace){ //Movie, movie_length, size_x, size_y, width, start, skip, traces, traces_out
	int cur_frame, cur_x, cur_y, j, wk, I;
	
	intensity_trace->clear();
	
		j = 0;		//frameindex in trace.... start with first 
		cur_frame = start; //frameindex in the movie
		cur_x = trace->at(j).at(1);
		cur_y = trace->at(j).at(2);
		while(cur_frame < movie_length){
			
			
			//intgerate here
			I = 0;
			for(int k=max(0,cur_x-w); k<=min(size_x-1,cur_x+w);k++){ 	
				wk = (int) (  sqrt(pow(w,2)-pow(cur_x-k,2))  ); //try circular
				for(int l=max(0,cur_y-wk); l<=min(size_y-1,cur_y+wk);l++){
					I += (int) (movie[size_x*size_y*cur_frame+size_y*k+l]);
					//mexPrintf("k=%i\tl=%i\twk=%i\n", k, l, I);
				}
			} 
			//mexPrintf("x=%i\ty=%i\tI=%i\n", cur_x, cur_y, I);
			vector<int> addthis;
			addthis.push_back(cur_frame);
			addthis.push_back(cur_x);
			addthis.push_back(cur_y);
			addthis.push_back(I);
			
			for(int k = 3; k < trace->at(j).size(); k++){ //the other variables

				addthis.push_back(trace->at(j).at(k));
			}  
			intensity_trace->push_back(addthis);
			
			//increment cur_frame
			cur_frame += skip+1;

			if(j+1 < trace->size()-1){
			if(cur_frame == trace->at(j+1).at(0)){
				j++;
				cur_x = trace->at(j).at(1);
				cur_y = trace->at(j).at(2);
				
			}
			}
			

		}



}




void PeakProcessing::integrate_double(double * movie, int movie_length, int size_x, int size_y, int w, int start, int skip,  std::vector< std::vector<double> >  * trace,  std::vector< std::vector<double> >  * intensity_trace){ //Movie, movie_length, size_x, size_y, width, start, skip, traces, traces_out
	int cur_frame, cur_x, cur_y, j, wk;
	double I = 0-0;
	
	intensity_trace->clear();
	
		j = 0;		//frameindex in trace.... start with first 
		cur_frame = start; //frameindex in the movie
		cur_x = (int) trace->at(j).at(1);
		cur_y = (int) trace->at(j).at(2);
		while(cur_frame < movie_length){
			
			
			//intgerate here
			I = 0;
			for(int k=max(0,cur_x-w); k<=min(size_x-1,cur_x+w);k++){ 	
				wk = (int) (  sqrt(pow(w,2)-pow(cur_x-k,2))  ); //try circular
				for(int l=max(0,cur_y-wk); l<=min(size_y-1,cur_y+wk);l++){
					I +=  (movie[size_x*size_y*cur_frame+size_y*k+l]);
					//mexPrintf("k=%i\tl=%i\twk=%i\n", k, l, I);
				}
			} 
			//mexPrintf("x=%i\ty=%i\tI=%i\n", cur_x, cur_y, I);
			vector<double> addthis;
			addthis.push_back(cur_frame);
			addthis.push_back(cur_x);
			addthis.push_back(cur_y);
			addthis.push_back(I);
			
			for(int k = 3; k < trace->at(j).size(); k++){ //the other variables

				addthis.push_back(trace->at(j).at(k));
			}  
			intensity_trace->push_back(addthis);
			
			//increment cur_frame
			cur_frame += skip+1;

			if(j+1 < trace->size()-1){
			if(cur_frame == (int) trace->at(j+1).at(0)){
				j++;
				cur_x = (int) trace->at(j).at(1);
				cur_y = (int) trace->at(j).at(2);
				
			}
			}
			

		}



}


void PeakProcessing::append_traces(double * movie,int movie_length,int size_x, int size_y, std::vector< std::vector< std::vector<int> > > * traces, std::vector< std::vector< std::vector<int> > > * itraces, double r_find, double h_min, double r_integrate, std::vector<int> * framenumber){
	int cur_frame, trace_index, w, wk, intensity;
	double x, y, d_min;
	vector< vector<int> > cur_peaks;
	set<int> not_assigned;
	set<int>::iterator it;

	w = r_find;
	//make sure framenumber = movie_length
	if (framenumber->size() != movie_length) {
		mexPrintf("ERROR: Size of framnumber-vector does not equal length of movie!\n");
		mexPrintf("%i framenumber \t %i movie_length", framenumber->size(), movie_length);
	}
	else {
		
		
		
	//init the raster
	int raster[size_x][size_y];

	for (int i = 0; i < size_x; i++) {
		for (int j = 0 ; j < size_y; j++) {
			raster[i][j] = -1;
		}
	}	
		
	// fill the raster
	for (int i=0; i<traces->size(); i++) {
		if (raster[ traces->at(i).back().at(1) ][ traces->at(i).back().at(2) ] != i && raster[ traces->at(i).back().at(1) ][ traces->at(i).back().at(2) ] != -1) {
			mexPrintf("Warning: trace is assigned to a raster spot, which is already assigned!\n");
		}
		for (int j=0; j < traces->at(i).size(); j++) {
			raster[ traces->at(i).at(j).at(1) ][ traces->at(i).at(j).at(2) ] = i;
			//mexPrintf("set %i,%i to %i\n", traces->at(i).at(j).at(1) , traces->at(i).at(j).at(2), i);
		}
	}
		
	
	

		//loop through movie, identify peaks and add the to traces and itraces
	for (int i=0; i < movie_length; i++) {
		not_assigned.clear();
		//fill not_assigend set
		for (int j=0; j<traces->size(); j++) {
			not_assigned.insert(j);
		}
		
		cur_peaks.clear();
		PeakProcessing::find_peaks(movie, size_x*size_y*i, size_x, size_y , r_find, h_min, &cur_peaks, 0);	//find peaks in current frame
		cur_frame = framenumber->at(i);
		//append peaks to traces
		for (int j=0; j <cur_peaks.size(); j++) {
			x = cur_peaks[j][0]; 
			y = cur_peaks[j][1];
			trace_index = -1;
			d_min = w+1;
			for(int k=max(0,x-w); k<=min(size_x-1,x+w);k++){ 	
				wk = (int) (  sqrt(pow(w,2)-pow(x-k,2))  ); //try circular
				for(int l=max(0,y-wk); l<=min(size_y-1,y+wk);l++){
					
					if ( (raster[k][l] != -1) && (sqrt(pow(k-x,2)+pow(l-y,2)) < d_min) ) {
						trace_index = raster[k][l];
						d_min = sqrt(pow(k-x,2)+pow(l-y,2));
					}
				}
			} 
			
			//integrate
			intensity = 0;
			for(int k=max(0,x-r_integrate); k<=min(size_x-1,x+r_integrate);k++){ 	
				wk = (int) (  sqrt(pow(r_integrate,2)-pow(x-k,2))  ); //try circular
				for(int l=max(0,y-wk); l<=min(size_y-1,y+wk);l++){
					intensity +=  movie[size_x*size_y*i+size_y*k+l];
				}
			} 
			
			if (trace_index >= 0) { //match
				//mexPrintf("match\n");
				//append to trace
				vector<int> tmp;
				tmp.push_back(cur_frame);
				tmp.push_back(x);
				tmp.push_back(y);
				traces->at(trace_index).push_back(tmp);
				//update raster 
				raster[(int)x][(int)y] = trace_index;
				//append to itraces
				vector<int> tmp2;
				tmp2.push_back(cur_frame);
				tmp2.push_back((int)x);
				tmp2.push_back((int)y);				
				tmp2.push_back(intensity);
				itraces->at(trace_index).push_back(tmp2);
				
				//remove from not_assigned list
				not_assigned.erase(trace_index);
			}
			else { // no match found	
				//mexPrintf("no match\n");

				raster[(int)x][(int)y] = traces->size(); //this will be a new trace
				//append to traces
				vector< vector<int> > tmp;
				vector<int> tmp2;
				tmp2.push_back(cur_frame);
				tmp2.push_back(x); 
				tmp2.push_back(y);
				tmp.push_back(tmp2);
				traces->push_back(tmp);
				
				//append to itraces
				vector< vector<int> > tmp4;
				vector<int> tmp3;
				tmp3.push_back(cur_frame);
				tmp3.push_back(x); 
				tmp3.push_back(y);
				tmp3.push_back(intensity);

				
				tmp4.push_back(tmp3);
				itraces->push_back(tmp4);
				
			}

			
		}
		
		// add frames to not assigned traces
		for (it=not_assigned.begin(); it != not_assigned.end(); it++) {
			trace_index = *it;
			x = traces->at(trace_index).back().at(1);
			y = traces->at(trace_index).back().at(2);
			//mexPrintf("%i\n", trace_index);
			//integrate
			intensity = 0;
			for(int k=max(0,x-w); k<=min(size_x-1,x+w);k++){ 	
				wk = (int) (  sqrt(pow(w,2)-pow(x-k,2))  ); //try circular
				for(int l=max(0,y-wk); l<=min(size_y-1,y+wk);l++){
					intensity +=  movie[size_x*size_y*i+size_y*k+l];
				}
			} 
			vector<int> tmp;
			tmp.push_back(cur_frame);
			tmp.push_back((int)x);
			tmp.push_back((int)y);				
			tmp.push_back((int)intensity);
			itraces->at(trace_index).push_back(tmp);
		}
		
		
	}
	
	}
	
}



void PeakProcessing::append_traces_to_position(double * movie,int movie_length,int size_x, int size_y, std::vector< std::vector< std::vector<int> > > * itraces, std::vector< std::vector<int> > * positions, double r_integrate, std::vector<int> * framenumber){
	int cur_frame, x, y;
	double I, wk, w;
	w = r_integrate;
	
	if (framenumber->size() != movie_length) {
		mexPrintf("ERROR: Size of framnumber-vector does not equal length of movie!\n");
		mexPrintf("%i framenumber \t %i movie_length", framenumber->size(), movie_length);
	}
	else {
		
		
		
		
	
	for (int i=0; i<movie_length; i++) {
		cur_frame = framenumber->at(i); 
		
		//loop over positions
		for (int j=0; j < positions->size(); j++) {
			x = positions->at(j).at(0);
			y = positions->at(j).at(1);

			I = 0; //integrate
			for(int k=max(0,x-w); k<=min(size_x-1,x+w);k++){ 	
				wk = (int) (  sqrt(pow(w,2)-pow(x-k,2))  ); //try circular
				for(int l=max(0,y-wk); l<=min(size_y-1,y+wk);l++){
					I +=  movie[size_x*size_y*i+size_y*k+l];
				}
			} 
			
			//append to itraces
			
			vector<int> tmp2;
			tmp2.push_back(cur_frame);
			tmp2.push_back(x);
			tmp2.push_back(y);
			tmp2.push_back(I);
			
			if (itraces->size() < positions->size()) { // if itraces is empty or not filled jet
				vector< vector<int> > tmp;
				tmp.push_back(tmp2);
				itraces->push_back(tmp);

			}
			else {
				itraces->at(j).push_back(tmp2); // if the trace already exists

			}

			
			
			
		}
		
	}
	
		
		
	}
}






void PeakProcessing::map_traces(std::vector< std::vector<double> > * pos1, std::vector< std::vector<double> > * pos2, std::vector< std::vector<double> > * pos3, std::vector< std::vector<int> > * map, double r){
	double x1,y1,x2,y2,x3,y3,d12,d13,d23, r_min;	
	int j_min, k_min, index2, index3;
    set<int> not_assigned2;
    set<int> not_assigned3;
    
    for(int i=0; i< pos2->size(); i++) { not_assigned2.insert(i);} //fill assigned vales
    for(int i=0; i< pos3->size(); i++) { not_assigned3.insert(i);}
    set<int>::iterator it2;
	set<int>::iterator it3;

	for (int i=0; i < pos1->size(); i++) {
		j_min = -1;
		k_min = -1;
		r_min = 3*r;
		x1 = pos1->at(i).at(0);
		y1 = pos1->at(i).at(1);
		for (it2 = not_assigned2.begin();it2 != not_assigned2.end(); it2++) {
            index2 = *it2;
			x2 = pos2->at(index2).at(0);
			y2 = pos2->at(index2).at(1);
			d12 = sqrt(pow(x1-x2,2) + pow(y1-y2,2));
			if (d12 < r) {
				for ( it3 = not_assigned3.begin(); it3 != not_assigned3.end(); it3++) {
                    index3 = *it3;
					x3 = pos3->at(index3).at(0);
					y3 = pos3->at(index3).at(1);
					d13 = sqrt(pow(x1-x3,2) + pow(y1-y3,2));
					d23 = sqrt(pow(x2-x3,2) + pow(y2-y3,2));
					//mexPrintf("d12 = %.4f\t d13 = %.4f \t d23 = %.4f\n", d12, d13, d23);

					if (d23 < r && d13 < r && (d12+d13+d23) < r_min){
						r_min = d12+d13+d23;
						j_min = index2;
						k_min = index3;	
						//mexPrintf("hit: %f %f\t%f %f\t%f %f \n", x1,y1,x2,y2,x3,y3);
						//mexPrintf("%i %i %i\n", i, j_min, k_min);
					}					
				}
			}
		}
		
		if (j_min >= 0 && k_min >= 0) {
			//mexPrintf("Added: %i %i %i\n", i, j_min, k_min);

			vector<int> tmp;
			tmp.push_back(i);
			tmp.push_back(j_min);
			tmp.push_back(k_min);
			map->push_back(tmp);
            
            //remove j_min and k_min from not_assigned2 and not_assigend_3
            not_assigned2.erase(j_min);
            not_assigned3.erase(k_min);
		}	
		
		
	}
	
	
	
	
}











