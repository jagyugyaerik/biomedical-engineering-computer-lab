%LABOR3 - alapfrekvencia eloszlása öt kiválasztott dallamban
close all

figure(1)
hist([F01, F02, F03, F04, F05],300:50:800)
xlabel('Alapfrekvencia értékek (Hz)')
ylabel('Elõfordulási gyakoriság')
title('Az öt kiválasztott dallam alapfrekvencia-eloszlása')

saveas(figure(1), 'lab3_7_hisztogram', 'jpg')

