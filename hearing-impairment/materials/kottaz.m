%
%KOTTÁZ függvény
%
% A beadott dallamot az ötvonalas ábrázolásmódban jeleníti meg.
%
% KOTTAZ(DALLAM): a DALLAM egy olyan vektor, mely a dallamot tartalmazza.
%
% KOTTAZ(DALLAM,IDOTENGELY): a beadott IDOTENGELY vektor alapján ábrázol.
%
% KOTTAZ(DALLAM,IDOTENGELY,MAX_IDO): a MAX_IDO-vel a maximális ábrázolási
% idõ adható meg, függetlenül a dallam hosszától
%
% KOTTAZ(DALLAM,MOD): ahol a MOD-ban megadhatjuk, hogy miként ábrázolja a
% dallamot
%       MOD = 1     az idoskála a pontok sorszámozásával egyezik meg (alap)
%       MOD = 2     az idoskála normalizált 0 és 1 között
%
%         
% Várallyay György - 2006.07.01.

function kottaz(dallam,idotengely_vagy_mod,max_ido)


% ----- bemeneti adatok ellenõrzése -----

switch nargin
    case 0
        error('HIBA! Nincs bemeneti adat!')
    case 1
        idotengely = 1:length(dallam);
    case 2
        if length(idotengely_vagy_mod) == 1
            % a MOD üzemmódban vagyunk
            mod = idotengely_vagy_mod;
            switch mod
                case 1
                    idotengely = 1:length(dallam);
                case 2
                    idotengely = linspace(0,1,length(dallam));
            end
        else
            % IDOTENGELY üzemmod
            idotengely = idotengely_vagy_mod;
        end
    case 3
        idotengely = idotengely_vagy_mod;
end

if not(size(dallam)==size(idotengely))
    error('HIBA! A DALLAM és az IDOTENGELY mérete nem azonos!')
end


% ----- ábrázolás -----

abra = semilogy(idotengely,dallam,'o');

set(abra,...
    'LineWidth',2,...
    'MarkerSize',8,...
    'MarkerFaceColor','b',...
    'MarkerEdgeColor','none'...
    )

v = [329.63 397.70 479.82 578.91 698.46];  % öt vonal pontos értéke
vh = [330 400 480 580 700];                % öt vonal hozzávetõleges értéke

set(gca,...
    'YTick',v,...
    'YTickMode','manual',...
    'YTickLabel','330|400|480|580|700',...
    'GridLineStyle','-',...
    'YGrid','on',...
    'YLim',[200 1000]...                   % a teljes megjelenített frek-
    )                                      % venciatartomány: 200-1000 Hz

xlabel('idötengely')
ylabel('frekvencia (Hz)')
title('Beadott dallam')

if nargin == 3
    set(gca,'Xlim',[0 max_ido])
end
