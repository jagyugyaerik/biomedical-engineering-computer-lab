function PVOL=readRawR32(filename, dim)
%----------------------------------
fid = fopen(filename,'rb');
if (size(dim)==[1 1])
	dim=[dim dim dim];
end;
if (fid < 0)
    PVOL=0;
else
    PVOL = fread(fid, prod(dim), 'single');
    PVOL = reshape(PVOL, dim);
    fclose(fid);
end;
