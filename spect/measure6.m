clear all;
close all;
%% Radon
% Feladat 1
img = zeros([128 128]);
img(60, 50) = 100;
P180 = radon(img, 0:3:177);
P360 = radon(img, 0:3:357);
showImage(img, 'Mérési kép', 1);
saveas(figure(1), 'results/meresikep', 'jpg');
showImage(P180, 'Mérési kép radon transzformáltja', 2);
saveas(figure(2), 'results/meresikepRadon180', 'jpg');
showImage(P360, 'Mérési kép radon transzformáltja', 3);
saveas(figure(3), 'results/meresikepRadon360', 'jpg');

% Feladat 2
close all;
for numberOfAngles = 1:20
    angleDiff = 180/numberOfAngles;
    proj = radon(img, 0:angleDiff:180);
    R = iradon(proj, 0:angleDiff:180, 'linear', 'none');
    showImage(R, ['Radon transzformáció rekonstrukció vetületek ', num2str(angleDiff), ' fokonként'], numberOfAngles);
    angle = replace(num2str(angleDiff), '.', ',');
    saveas(figure(numberOfAngles), ['results/iradon/iradonNumOfAngles',angle], 'jpg');
end

%% SM
% Feladat 1
img = zeros([64 64]);
img(21, 30) = 100;
A = SPECT2DSystemMatrix(64, 0:3:359, 0);
PSM = SPECT2DForwardProjSM(img, A);
showImage(PSM, 'Radon trasnszformáció rendszermátrixszal', 1);
saveas(figure(1), 'results/sm/redszerMatrix', 'jpg');
Aerror = SPECT2DRRSystemMatrix(64, 0:3:359, 0, 0.97030, 0.017239, 2.0, 280, 6.0);
Perror = SPECT2DForwardProjSM(img, Aerror);
showImage(Perror, 'Radon transzformáció távolságfüggő elkenéses rendszermátrixszal', 2);
saveas(figure(2), 'results/sm/redszerMatrixTavolsagelkenes', 'jpg');
projRadon = radon(img, A.theta);
showImage(projRadon, 'Radon transzformáció', 3);
saveas(figure(3), 'results/sm/radon', 'jpg');

%% MLEM
close all
% Feladat 1
img = readRawR32('phantom-shepplogan-64x64.raw', [64 64]);
showImage(img, 'Shepp-Logan fantom', 1);
saveas(figure(1), 'results/mlem/mlemRaw', 'jpg');

S = SPECT2DSystem(64, 0:3:359, 0, 0.97030, 0.017239, 2.0, 280, 6.0);
proj = SPECT2DForwardProj(img, S);
showImage(proj, 'Shepp-Logan előrevetített kép', 2);
saveas(figure(2), 'results/mlem/mlemForward', 'jpg');

projRadon = SPECTRadon(img, 0:3:359);
showImage(projRadon, 'Radon transzformáció', 3);
saveas(figure(3), 'results/mlem/mlemRadon', 'jpg');

dL2MLEM = zeros(11,1);
dL2OSEM = zeros(11,1);
dCCMLEM = zeros(11,1);
dCCOSEM = zeros(11,1);
%%
rMLEM = SPECT2DMLEM(proj, S, 10);
dCCMLEM(1) = ImageDistance_CC( img, rMLEM );
dL2MLEM(1) = ImageDistance_L2(img, rMLEM);

rOSEM = SPECT2DOSEM(proj, S, 10, 4);
dL2OSEM(1) = ImageDistance_L2( img, rOSEM );
dCCOSEM(1) = ImageDistance_L2( img, rOSEM );

showImage(rMLEM, 'MLEM rekonsturkció iteráció száma: 10', 10, 122);
showImage(rOSEM, 'OSEM rekonsturkció iteráció száma: 10 subset: 4', 10, 121);
saveas(figure(10), 'results/mlem/osem10', 'jpg')
%%
for iter = 11:20
    rMLEM = SPECT2DMLEM(proj, S, iter, rMLEM);
    dCCMLEM(iter-10+1) = ImageDistance_CC( img, rMLEM );
    dL2MLEM(iter-10+1) = ImageDistance_L2( img, rMLEM );

    rOSEM = SPECT2DOSEM(proj, S, 10, 4, rOSEM);
    dL2OSEM(iter-10+1) = ImageDistance_L2( img, rOSEM );
    dCCOSEM(iter-10+1) = ImageDistance_CC(img, rOSEM);
    
    showImage(rMLEM, ['MLEM rekonsturkció iteráció száma: ', num2str(iter)], iter, 121);
    showImage(rOSEM, ['OSEM rekonsturkció iteráció száma: ', num2str(iter), ' subset: 4'], iter, 122);
    saveas(figure(iter), ['results/mlem/osem', num2str(iter)], 'jpg')
end

%% ATT
close all;
source = readRawR32('raw/phantom-ring-128x128.raw', [128 128]);
showImage(source, 'Fantom gyűrű', 1);
mumap = readRawR32('raw/attenuation-cyl-128x128.raw', [128 128]);
S = SPECT2DSystem(128, 0:3:359, 128, 0.97030, 0.017239, 2, 280, 3 );
proj = SPECT2DAttenuatedForwardProj( source, S, mumap );
proj = SPECT2DAddNoise( proj );
rOSEM = SPECT2DMLEM( proj, S, 10, 4 );
showImage(rOSEM, 'OSEM rekonsturkció iteráció száma: 10 subset: 4', 2);
rIradon = iradon(proj, 0:3:359, 'linear', 'none');
showImage(rIradon, 'Iradon rekonstrukció', 3);