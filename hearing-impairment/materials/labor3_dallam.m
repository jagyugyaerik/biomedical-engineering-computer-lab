%LABOR3 - dallamk�sz�t�s
close all

display('**********************')
display('LABOR3 - dallamelemz�s')
display('**********************')

display('Adja meg a kiv�lasztott szegmens sorsz�m�t, valamint kezd�- �s v�gpontj�t!')
sorszam = str2num(input('A kiv�lasztott szegmens (eredeti) sorsz�ma: ','s'));
kezdet = str2num(input(['A(z) ',int2str(sorszam),'-ik szegmens kezd�pontja: '], 's'));
veg = str2num(input(['A(z) ',int2str(sorszam),'-ik szegmens v�gpontja: '], 's'));

yy = y(crySegment(1,1):crySegment(1,2));             % az aktu�lis szegmens kiv�g�sa

N = 4096;                       % ablakm�ret: ennyi mint�b�l �ll egy ablak
f = (0:N/2)*fs/N;               % a spektrum frekvenciatengely�nek l�trehoz�sa

Nwin = floor(length(yy)/N);     % ennyi ablak f�r bele a szegmensbe

clear F0
figure(1)
Nwin = floor(length(yy)/N);
for i=1:Nwin
	win = yy((i-1)*N+1:i*N);
	W = abs(fft(win,N));
	plot(f,W(1:N/2+1))
	set(gca,'XLim',[0 4000])
	title(['Ablak sorsz�ma: ',int2str(i)])
	[F0_i,I_i] = ginput(1);
	F0(i) = F0_i;
end

T = (1:Nwin)*(Ts*N);

figure(2)
kottaz(F0,T)
title(['A ',siras_neve,' f�jl ',int2str(sorszam),'-ik szegmens�nek dallama'])

if exist('F04')
    F05 = F0;
    display('�t�dik kiv�lasztott szegmens elt�rolva F05-ben')
    saveas(figure(2), 'lab3_6_dallam', 'jpg')
elseif exist('F03')
    F04 = F0;
    display('Negyedik kiv�lasztott szegmens elt�rolva F04-ben')
    saveas(figure(2), 'lab3_5_dallam', 'jpg')
elseif exist('F02')
    F03 = F0;
    display('Harmadik kiv�lasztott szegmens elt�rolva F03-ben')
    saveas(figure(2), 'lab3_4_dallam', 'jpg')
elseif exist('F01')
    F02 = F0;
    display('M�sodik kiv�lasztott szegmens elt�rolva F02-ben')
    saveas(figure(2), 'lab3_3_dallam', 'jpg')
else
    F01 = F0;
    display('Els� kiv�lasztott szegmens elt�rolva F01-ben')
    saveas(figure(2), 'lab3_2_dallam', 'jpg')
end
    
