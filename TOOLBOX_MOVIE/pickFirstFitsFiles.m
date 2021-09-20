function [ files ] = pickFirstFitsFiles( pname, channel )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    tmp = dir([pname filesep channel '*.fits']);
    tmp = {tmp.name}; % make a cell out of tmp
    
    expression = '_X[0-9]+.fits'; % expression to sort out X
    
    index_of_exp = regexp(tmp, expression);
    
    files = cell(1, 0);
    for i=1:length(tmp)
        if isempty(index_of_exp{i})
            files = [files tmp{i}];
        end
    end
    
    
end

