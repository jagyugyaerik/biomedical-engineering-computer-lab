%LABOR3 - �tlagos szegmenshossz �s szegmensgyakoris�g kisz�m�t�sa
% meg kell adni a szegmensek kezd�- �s v�gpontjait
% ha nincs t�bb szegmens, �rj be 0-�t mindk�t helyre!
clear start_end

display('****************************************************************')
display('LABOR3 - �tlagos szegmenshossz �s szegmensgyakoris�g kisz�m�t�sa')
display('****************************************************************')

kezdet = '999';
num_seg = 0;

while str2double(kezdet)~=0
    num_seg = num_seg+1;
    kezdet = input(['A(z) ',int2str(num_seg),'-ik szegmens kezdete: '], 's');
    veg = input(['A(z) ',int2str(num_seg),'-ik szegmens v�ge: '], 's');
    start_end(num_seg,:) = [str2double(kezdet),str2double(veg)];
end
start_end(end,:) = [];

display(['�sszesen:             ',int2str(num_seg-1),' szegmens'])

Ttot = length(y)*Ts;
Tseg = sum((start_end(:,2)-start_end(:,1))*Ts);

display(['Szegmensgyakoris�g =  ',num2str(0.1*round(1000*Tseg/Ttot)),' %'])
display(['�tlagos szegmensid� = ',num2str(mean((start_end(:,2)-start_end(:,1))*Ts)), 'sec'])
