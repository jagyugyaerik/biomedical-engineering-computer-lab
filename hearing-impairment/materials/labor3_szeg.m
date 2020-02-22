%LABOR3 - szegmens(ek) meghat�roz�sa

display('***********************************')
display('LABOR3 - szegmens(ek) meghat�roz�sa')
display('***********************************')

figure(2)
plot(y)                 % s�r�s �br�zol�sa a mint�k sorsz�m�val
title(fileName)
xlabel('minta sorszama'), ylabel('amplitudo')

[x1,y1] = ginput(1);		% kattintsunk a szegmens elej�re
[x2,y2] = ginput(1);		% kattintsunk a szegmens v�g�re

display(['A szegmens kezd�pontja: ',int2str(x1)])
display(['A szegmens v�gpontja:   ',int2str(x2)])
display(['A szegmens id�tartama:  ',num2str((x2-x1+1)*Ts), ' sec'])