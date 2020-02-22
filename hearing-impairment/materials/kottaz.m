%
%KOTT�Z f�ggv�ny
%
% A beadott dallamot az �tvonalas �br�zol�sm�dban jelen�ti meg.
%
% KOTTAZ(DALLAM): a DALLAM egy olyan vektor, mely a dallamot tartalmazza.
%
% KOTTAZ(DALLAM,IDOTENGELY): a beadott IDOTENGELY vektor alapj�n �br�zol.
%
% KOTTAZ(DALLAM,IDOTENGELY,MAX_IDO): a MAX_IDO-vel a maxim�lis �br�zol�si
% id� adhat� meg, f�ggetlen�l a dallam hossz�t�l
%
% KOTTAZ(DALLAM,MOD): ahol a MOD-ban megadhatjuk, hogy mik�nt �br�zolja a
% dallamot
%       MOD = 1     az idosk�la a pontok sorsz�moz�s�val egyezik meg (alap)
%       MOD = 2     az idosk�la normaliz�lt 0 �s 1 k�z�tt
%
%         
% V�rallyay Gy�rgy - 2006.07.01.

function kottaz(dallam,idotengely_vagy_mod,max_ido)


% ----- bemeneti adatok ellen�rz�se -----

switch nargin
    case 0
        error('HIBA! Nincs bemeneti adat!')
    case 1
        idotengely = 1:length(dallam);
    case 2
        if length(idotengely_vagy_mod) == 1
            % a MOD �zemm�dban vagyunk
            mod = idotengely_vagy_mod;
            switch mod
                case 1
                    idotengely = 1:length(dallam);
                case 2
                    idotengely = linspace(0,1,length(dallam));
            end
        else
            % IDOTENGELY �zemmod
            idotengely = idotengely_vagy_mod;
        end
    case 3
        idotengely = idotengely_vagy_mod;
end

if not(size(dallam)==size(idotengely))
    error('HIBA! A DALLAM �s az IDOTENGELY m�rete nem azonos!')
end


% ----- �br�zol�s -----

abra = semilogy(idotengely,dallam,'o');

set(abra,...
    'LineWidth',2,...
    'MarkerSize',8,...
    'MarkerFaceColor','b',...
    'MarkerEdgeColor','none'...
    )

v = [329.63 397.70 479.82 578.91 698.46];  % �t vonal pontos �rt�ke
vh = [330 400 480 580 700];                % �t vonal hozz�vet�leges �rt�ke

set(gca,...
    'YTick',v,...
    'YTickMode','manual',...
    'YTickLabel','330|400|480|580|700',...
    'GridLineStyle','-',...
    'YGrid','on',...
    'YLim',[200 1000]...                   % a teljes megjelen�tett frek-
    )                                      % venciatartom�ny: 200-1000 Hz

xlabel('id�tengely')
ylabel('frekvencia (Hz)')
title('Beadott dallam')

if nargin == 3
    set(gca,'Xlim',[0 max_ido])
end
