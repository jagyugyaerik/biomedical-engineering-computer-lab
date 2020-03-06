function saveRawR32(filename, vol)
%----------------------------------
fid = fopen(filename,'wb');
if (fid < 0)
    printf('Couldn''t open file for writing.\n');
else
    fwrite(fid, vol, 'single');
    fclose(fid);
end;
