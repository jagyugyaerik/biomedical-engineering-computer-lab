%LABOR3 - inicializ�l�s, s�r�sf�jl beolvas�sa, �br�zol�sa

close all
clear all
clc

display('********************************************************')
display('LABOR3 - inicializ�l�s, s�r�sf�jl beolvas�sa, �br�zol�sa')
display('********************************************************')

siras_neve = input('S�r�sf�jl neve: ','s')

[y,fs] = audioread(siras_neve);
y = y-mean(y);          % DC komponens elt�vol�t�sa
Ts = 1/fs;              % mintav�teli id�(k�z)

display(['Az elemzend� s�r�sf�jl neve:  ',siras_neve])
display(['Mintav�teli frekvencia (fs) = ',int2str(fs),' Hz'])
display(['A s�r�s pontjainak sz�ma:     ',int2str(length(y))])
display(['A s�r�s hossza (Ttot) =       ',num2str(length(y)*Ts),' sec'])

display('FIGURE1: a s�r�s �br�zol�sa a mint�k sorsz�m�val')
figure(1)
plot(y)                 % s�r�s �br�zol�sa a mint�k sorsz�m�val
title(siras_neve)
xlabel('minta sorszama'), ylabel('amplitudo')


t = (1:length(y))/fs;	% id�tengely l�trehoz�sa

display('FIGURE2: a s�r�s �br�zol�sa id�tengellyel')
figure(2)
plot(t,y)               % s�r�s �br�zol�sa id�tengellyel
title(siras_neve)
xlabel('ido [s]'), ylabel('amplitudo')

saveas(figure(2), 'lab3_1_siras_idojele', 'jpg')
