%LABOR3 - alapfrekvencia eloszl�sa �t kiv�lasztott dallamban
close all

figure(1)
hist([F01, F02, F03, F04, F05],300:50:800)
xlabel('Alapfrekvencia �rt�kek (Hz)')
ylabel('El�fordul�si gyakoris�g')
title('Az �t kiv�lasztott dallam alapfrekvencia-eloszl�sa')

saveas(figure(1), 'lab3_7_hisztogram', 'jpg')

