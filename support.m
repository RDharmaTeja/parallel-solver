g = load('Rchannel-3D-2.txt');
x = 181;y = 65; z = 2;
g = g(2:end,:);
g1 = [];
g2 = [];
for i = 0:y*z-1
g1 = [g1; g(i*x +1: i*x+91,:)];
g2 = [g2; g(i*x +92: i*x+181,:)];
end
g1 = [91 65 2; g1];
g2 = [90 65 2; g2];
save('process_1.txt', 'g1', '-ASCII','-append');
save('process_2.txt', 'g2', '-ASCII','-append');