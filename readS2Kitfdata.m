function f = readS2Kitfdata(fname)
% usage
%        f = readS2Kitfdata(fname)
% reads in the interleaved function
% data used by S2Kit, which is in raster scanned format
% obtains BW from the 
% the output f is a 2BW x 2BW matrix


if ( nargin < 1 )
    [fname,pname]=uigetfile('*.dat','Pick a data file');
    fname = [pname,fname];
end;
    
fin = fopen(fname,'rt');
% read in
Fraw = fscanf(fin,'%f');
fclose(fin);
% estimate bw
BW = sqrt(length(Fraw)/8); 
% if complex data, we have 2*4B^2 = 8B^2 data interleaved
% hence 8B^2 = L ==> B = sqrt(L/8);
% if real data we have 4B^2
% hence B = sqrt(L/4);
if (BW ~= fix(BW)) % not an int, must be real data
    BW = sqrt(length(Fraw)/4);
    Fcomplex = Fraw;  
else % otherwise it was interleaved, so del-interleavit
    Fcomplex = Fraw(1:2:end)+j*Fraw(2:2:end);
end;
f = reshape(Fcomplex,2*BW,2*BW);
% matlab uses columnwise, put into row-wise
f = conj(f');