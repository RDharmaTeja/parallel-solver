g = load('main.txt');
x = 181;y = 65; z = 2;
cut = 100;
g = g(2:end,:);
g1 = [];
g2 = [];
for i = 0:y*z-1
g1 = [g1; g(i*x +1: i*x+cut,:)];
g2 = [g2; g(i*x +cut: i*x+181,:)];
end
g1 = [cut 65 2; g1];
g2 = [181-cut+1 65 2; g2];
save('process_0.txt', 'g1', '-ASCII','-append');
save('process_1.txt', 'g2', '-ASCII','-append');