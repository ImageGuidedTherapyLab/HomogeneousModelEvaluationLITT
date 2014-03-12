% This function is shamelessly copied from Chris MacLellan's
% /FUS4/data2/CJM/matlab/readVTK.m   but has a few edits to make its use
% lucid, including an error message if fopen doens't work.

% E.g. input   array = readVTK_SJF('temperature',161);
% For an series beginning with temperature.0000.vtk and ending with 0161.

function array = readVTK_SJF(vtkfile,timepoints);
%for VTK files starting with 0000. e.g.
%VTKfile.0000.vtk,VTKfile.0001.vtk...
%NOT
%VTKfile.0001.vtk,VTKfile.0002.vtk....

array=zeros(256);
for ii = 1:timepoints


%fid = fopen(vtkfile,'r','b');
fid = fopen(sprintf('%s.%04d.vtk',vtkfile,ii-1),'r','b');

% if fid == -1
%     disp('Error, could not open file')
%     return
% end

fgetl(fid); % # vtk Datafile Version 3.6
fgetl(fid); % comments
fgetl(fid); % BINARY
fgetl(fid); % DATASET STRUCTURED_POINTS

s = fgetl(fid); % DIMENSIONS nx ny nz
sz = sscanf(s, '%*s%d%d%d').';

fgetl(fid); % SPACING
fgetl(fid); % ORIGIN
fgetl(fid); %POINT_DATA

fgetl(fid); % SCALARS name data_type (eg. SCALARS scalars float)
%svstr = sscanf(s, '%s', 1)
%dtstr = sscanf(s, '%*s%*s%*s')
fgetl(fid); % the LOOKUP_TABLE

V = fread(fid,prod(sz),'float');

V = reshape(V,sz);
array(:,:,ii) = V(:,:);

fclose(fid);
end
