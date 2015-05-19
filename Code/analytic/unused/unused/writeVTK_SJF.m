function [] = writeVTK_SJF(vol,vtkfile, header)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Usage: writeVTK(vol,vtkfile)
%
%   vol:     The 3D matrix to be saved to file
%   vtkfile: The output filename (string)
%   notes:   Only writes binary STRUCTURED_POINTS
%  
% Erik Vidholm 2006

% edited by Josh Yung 5/4/10

% edited by SJ Fahrenholtz 12/3/2013
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ischar(vtkfile) == 0
    error('VTKFILE must be a string.')
end

% dimensions
volinfo = whos('vol');

sz = volinfo.size;

X = sz(1); Y = sz(2); Z = 1;
if( length(sz) == 4 )
      Z = sz(4); 
end
if length(sz) == 2
    sz(3) = 1;
end


% pixel spacing
spX = header.PixelSpacing(1);
spY = header.PixelSpacing(2);
spZ = header.SliceThickness;

% origin
oX = header.ImagePositionPatient(1);
oY = header.ImagePositionPatient(2);
oZ = header.ImagePositionPatient(3);

for ii = 1:sz(3)
    if length(sz) == 4
        voltemp = squeeze(vol(:,:,ii,:));
    else
        voltemp = squeeze(vol(:,:,ii));
    end

%vtkfile = 'vtk
fid = fopen(sprintf('%s.%04d.vtk',vtkfile,ii-1),'wt','b');

% write header
fprintf(fid, '%s\n', '# vtk DataFile Version 3.6');
fprintf(fid, '%s\n', 'created by writeVTK (Matlab implementation by Erik Vidholm)');
fprintf(fid, '%s\n', 'BINARY');  
fprintf(fid, '%s\n', 'DATASET STRUCTURED_POINTS');  
fprintf(fid, '%s%d%c%d%c%d\n', 'DIMENSIONS ', X, ' ', Y, ' ', Z);
fprintf(fid, '%s%.6f%c%.6f%c%d\n', 'SPACING ', spX, ' ', spY, ' ', spZ); 
fprintf(fid, '%s%f%c%f%c%f\n', 'ORIGIN ', oX, ' ', oY, ' ', oZ); 
fprintf(fid, '%s%d\n', 'POINT_DATA ', X*Y*Z);
fprintf(fid, 'SCALARS scalars float 1\n');

% tp = volinfo.class;
% if( strcmp(tp, 'uint8') > 0 )
%   fprintf(fid, '%s\n', 'SCALARS image_data unsigned_char');
% elseif( strcmp(tp, 'int8') > 0 )
%   fprintf(fid, '%s\n', 'SCALARS image_data char');
% elseif( strcmp(tp, 'uint16') > 0 )
%   fprintf(fid, '%s\n', 'SCALARS image_data unsigned_short');
% elseif( strcmp(tp, 'int16') > 0 )
%   fprintf(fid, '%s\n', 'SCALARS image_data short');
% elseif( strcmp(tp, 'uint32') > 0 )
%   fprintf(fid, '%s\n', 'SCALARS image_data unsigned_int');
% elseif( strcmp(tp, 'int32') > 0 )
%   fprintf(fid, '%s\n', 'SCALARS image_data int');
% elseif( strcmp(tp, 'single') > 0 )
%   fprintf(fid, '%s\n', 'SCALARS image_data float');
% elseif( strcmp(tp, 'double') > 0 )
%   fprintf(fid, '%s\n', 'SCALARS image_data float 1');
% end

fprintf(fid, '%s\n', 'LOOKUP_TABLE default');

% write data as binary
fwrite(fid,voltemp,'float');

% close file
fclose(fid);
end

