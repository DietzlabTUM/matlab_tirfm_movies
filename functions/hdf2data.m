% convert hdf5 file into data matrix

f = msgbox('Choose the folder with .hdf5 file');
uiwait(f);

pname=uigetdir(data_dir,'Choose the folder with .hdf5 file.');
cd (pname)

raw=hdf2mat;

path_out = [pname filesep datestr(now, 'yyyy-mm-dd_HH-MM') '_analysis'];
mkdir(path_out)

frames = max(raw(:,1));
spots = raw(end,end)+1;

data_tmp = cell(spots,1);
data = cell(spots,1);
for i=1:spots
    data_tmp{i}=zeros(frames,2);
end

for i = 1:spots
    
 tmp = raw(raw(:,end) == i-1, :,:);
     for j = 1:frames
         if tmp(tmp(:,1) == j)
             lines = tmp(tmp(:,1) == j, :,:);
             maxphoton = max(lines(:,4));
             data_tmp{i}(j,:) = lines(lines(:,4) == maxphoton , 2:3);
         else
             data_tmp{i}(j,:) = NaN;
         end
     end
    data{i}.fit=data_tmp{i};
         
end

cd (path_out)
save 'data.mat' data frames spots
