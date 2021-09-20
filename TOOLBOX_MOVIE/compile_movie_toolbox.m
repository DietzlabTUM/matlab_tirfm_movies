display('Compiling movie toolbox... please wait');
mex find_peaks2d.cpp peakprocessing.cpp;

mex append_traces.cpp peakprocessing.cpp;
mex append_traces_to_position.cpp peakprocessing.cpp;
mex map_traces.cpp peakprocessing.cpp;

%mex average_pos.cpp peakprocessing.cpp;
%mex make_traces.cpp peakprocessing.cpp;
%mex integrate_trace_double.cpp peakprocessing.cpp;
%mex integrate_traces.cpp peakprocessing.cpp;

%mex merge_2_traces.cpp peakprocessing.cpp;
%mex merge_traces.cpp peakprocessing.cpp;

display('Compiling done');



%{
%NLF functions
mex nlf.cpp peakprocessing.cpp;
mex nlf_alex.cpp peakprocessing.cpp;
mex nlf_fret.cpp peakprocessing.cpp;


%}