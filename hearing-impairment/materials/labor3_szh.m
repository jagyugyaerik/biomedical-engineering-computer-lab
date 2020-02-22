%LABOR3 - átlagos szegmenshossz és szegmensgyakoriság kiszámítása
% meg kell adni a szegmensek kezdõ- és végpontjait
% ha nincs több szegmens, írj be 0-át mindkét helyre!
clear start_end

display('****************************************************************')
display('LABOR3 - átlagos szegmenshossz és szegmensgyakoriság kiszámítása')
display('****************************************************************')

kezdet = '999';
num_seg = 0;

while str2double(kezdet)~=0
    num_seg = num_seg+1;
    kezdet = input(['A(z) ',int2str(num_seg),'-ik szegmens kezdete: '], 's');
    veg = input(['A(z) ',int2str(num_seg),'-ik szegmens vége: '], 's');
    start_end(num_seg,:) = [str2double(kezdet),str2double(veg)];
end
start_end(end,:) = [];

display(['Összesen:             ',int2str(num_seg-1),' szegmens'])

Ttot = length(y)*Ts;
Tseg = sum((start_end(:,2)-start_end(:,1))*Ts);

display(['Szegmensgyakoriság =  ',num2str(0.1*round(1000*Tseg/Ttot)),' %'])
display(['Átlagos szegmensidõ = ',num2str(mean((start_end(:,2)-start_end(:,1))*Ts)), 'sec'])
