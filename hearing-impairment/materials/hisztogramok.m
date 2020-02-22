% Hisztogramok
% Labor3 vizsgafeladat

h(1,:) = [9     0     2     7     3     0     0     0     0     0     0];
h(2,:) = [0     0     4    10    10     0     0     0     0     0     0];
h(3,:) = [0     0     0     0     0     0     0     0     0     0    21];
h(4,:) = [4    10    10     0     0     0     0     0     0     0     0];
h(5,:) = [0     0     0     0     0     1     4     8     9     2     3];
h(6,:) = [0     2    12    30    30     9     4     1     0     0     0];

skala = 300:50:800;

[r,ind] = sort(rand(6,1));

figure(99)
for i=1:6
    subplot(2,3,i)
    bar(skala,h(ind(i),:))
    set(gca,'XTick',300:100:800)
    grid on
end