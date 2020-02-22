clear all;
close all;
%% Init
% [file, path] = uigetfile('*.wav');
path = "/home/erik/code/Matlab/biomedical-engineering-computer-lab/hearing-impairment/raw/";
fileName = "CR0013.wav";

[y,fs] = audioread(fullfile(path, fileName));
y = y-mean(y);
Ts = 1/fs;

t = (1:length(y))/fs;

figure(1);
plot(y);
title(fileName);
xlabel('adat pont');
ylabel('amplitudo');
saveas(figure(1), "results/loadedFileAgainstDataPoints", "jpg");

figure(2);
plot(t, y);
title(fileName);
xlabel('idő [s]');
ylabel('amplitudo');
saveas(figure(2), "results/loadedFileAgainstTime", "jpg");
%% Extract Segments
close all;

figure(1);
plot(t, y);
title(fileName);
xlabel('idő [s]');
ylabel('amplitudo');
disp('Az ábra alapján add meg a szegmensek darabszámát!');
numberOfSegments = input("Szegmensek száma: ");

close all;
figure(1);
plot(y);
title(fileName);
xlabel('adat pont');
ylabel('amplitudo');

crySegmentsPositions = zeros(numberOfSegments, 2);
for segmentIndex = 1:numberOfSegments
    [x1,y1] = ginput(1);
    [x2,y2] = ginput(1);
    crySegmentsPositions(segmentIndex, 1) = x1;
    crySegmentsPositions(segmentIndex, 2) = x2;
end

% TODO refacotor to oneline
clear x1;
clear x2;
clear y1;
clear y2;


for segment=1:length(crySegmentsPositions)
   figure(segment);
   tseg = (1:length(y(crySegmentsPositions(segment,1):crySegmentsPositions(segment,2))))/fs;
   plot(tseg, y(crySegmentsPositions(segment,1):crySegmentsPositions(segment,2)));
   % title(fileName);
   title(strcat(fileName, " ", num2str(segment), ". szegmens"));
   xlabel('idő [s]');
   ylabel('amplitudo');
   saveas(figure(segment), strcat("results/segment", num2str(segment)), "jpg");
end

clear tseg;

%%  Calculations
tTotal = length(y) * Ts;
tSegmentsLength = zeros(1, length(crySegmentsPositions));
for segmentIndex = 1:length(tSegmentsLength)
    tSegmentsLength(segmentIndex) = (crySegmentsPositions(segmentIndex, 2) - crySegmentsPositions(segmentIndex, 1)) * Ts;
end

frequencySegment = sum(tSegmentsLength) / tTotal * 100;
meanSegment = sum(tSegmentsLength) / length(tSegmentsLength);


%% Frequency analysis
close all;

pointsOfWindow = 4096;
frequencyAxis = (0:pointsOfWindow/2)*fs/pointsOfWindow;

% for segmentIndex = 1:length(tSegmentsLength)
for segmentIndex = 1:3
    segment = y(crySegmentsPositions(segmentIndex, 1):crySegmentsPositions(segmentIndex, 2));
    numberOfWindows = floor(length(segment)/pointsOfWindow);
    baseFrequency = zeros(1, numberOfWindows);
    for windowIndex = 1:numberOfWindows
         window = segment((windowIndex-1) * pointsOfWindow + 1:windowIndex * pointsOfWindow);
         Y = abs(fft(window, pointsOfWindow));
         figure(1)
         plot(frequencyAxis, Y(1:pointsOfWindow/2+1));
         set(gca,'XLim',[0 4000]);
         title([int2str(windowIndex),' . ablak']);
         [x1, y1] = ginput(1);
         baseFrequency(windowIndex) = x1;
    end
    T = (1:numberOfWindows)*(Ts*pointsOfWindow);
    kottaz(baseFrequency,T)
    title([fileName ,int2str(segmentIndex),'. szegmensnek dallama'])
    disp('Nyomj entert a folytatéshoz...')
    pause
    
    m = 5;
    n = ceil(numberOfWindows / m);
    for windowIndex = 1:numberOfWindows
         figure(2)
         subplot(m,n,windowIndex)
         plot(frequencyAxis, Y(1:pointsOfWindow/2+1));
         title([int2str(windowIndex),' . ablak']);
         set(gca,'XLim',[0 4000]);
         xlabel('frekvencia [Hz]');
         ylabel('amplitudo');
    end
    saveas(figure(2), strcat("results/frequencyAnalysis", num2str(segmentIndex)), "jpg");
    close all;
end



% TODO refactor
clear f;
clear N;
clear x1;
clear y1;
