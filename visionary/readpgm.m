function image = readpgm(filename) ;

% READPGM  Reads a grayscale image file in pgm-format.
%          READPGM('filename.pgm') returns a matrix containing the 
%	   pixel values in the file PLUS ONE, since pgm-files
%          contains numbers 0..255 and Matlab works best with
%          numbers in the range 1..256.
%          
%          
%
%          See also SAVEPGM.

if nargin < 1, error('Too few input arguments.'); end
if ~isstr(filename), error('Argument must be a string.'); end


frid = fopen(filename,'r');
if frid == -1, error('Can not open file.'); end
idline = fgetl(frid) ;
if idline == -1, error('File is empty.'); end
if idline(1) ~= 'P', error('This is not a pbm, pgm or ppm-file.') ; end
if idline(2) ~= '5' & idline(2) ~= '2', error('This is not a pgm-file.') ; end

%disp('Warning, new version, see "help readpgm".');

whline = fgetl(frid);
if length(whline) == 0, whline(1) = '#'; end
while whline(1) == '#', 
  whline = fgetl(frid);
  if length(whline) == 0, whline(1) = '#'; end
end
[vek, count] = sscanf(whline,' %d',2);
width = vek(1);

%special fix by KAHL!
if count==1,
    whline = fgetl(frid);
    if length(whline) == 0, whline(1) = '#'; end
    while whline(1) == '#', 
      whline = fgetl(frid);
      if length(whline) == 0, whline(1) = '#'; end
    end
    [vek, count] = sscanf(whline,' %d',1);
    height = vek(1);
else %normal case
    height = vek(2);
end



whline = fgetl(frid);
if length(whline) == 0, whline(1) = '#'; end
while whline(1) == '#', 
  whline = fgetl(frid);
  if length(whline) == 0, whline(1) = '#'; end
end
[range, count] = sscanf(whline,' %d',1);

if idline(2) == '5',
  [image,count] = fread(frid,[width height],'uchar') ;
elseif idline(2) == '2',
  [image,count] = fscanf(frid,' %d ',[width height]);
end

image = image'+1; % Add one to avoid zeroes in Matlab
fclose(frid);
