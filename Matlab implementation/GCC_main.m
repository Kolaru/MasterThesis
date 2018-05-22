clear; clc;
%It will take 10 hours for number_of_nodes = 10000; times = 100;
number_of_nodes = 1000;
times = 10; 
N = 100; 
step = 0.02;
result = zeros(N,1);
for i = 1:N
    i
    for j = 1:times        
        c = step*i;
        A = ER(number_of_nodes,c);
        B = largestcomponent(A);
        result(i,j) = size(B,2);
    end
end
figure(1); hold on;
plot(step:step:step*N,mean(result/number_of_nodes,2),'*')