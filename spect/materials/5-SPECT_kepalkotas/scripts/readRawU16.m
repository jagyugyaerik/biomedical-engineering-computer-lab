function PVOL=readRawU16(filename, dim)
%----------------------------------
fid = fopen(filename,'r');
if (size(dim)==[1 1])
	dim=[dim dim dim];
end
if (fid < 0)
    PVOL=0;
else
    PVOL = fread(fid, prod(dim), 'uint16');
    if (prod(size(PVOL)) ~= prod(dim))
        PVOL=0;
    else
        PVOL = reshape(PVOL, dim);
    end;
    fclose(fid);
end;
