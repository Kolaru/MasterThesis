clear; clc;
load power.mat;
data = power;
max_id = max(max(data));

% Adjacency matrix
A = sparse(data(:,1),data(:,2),1);
A(max_id, max_id) = 0 ; % Make sure the matrix is a square by specifying the bottom right element
A = A + A';

%figure(1);
%h1 = plot(graph(A));
%layout(h1,'force','UseGravity',true)
%% shuffling
data_random = [data(:,1),data(randperm(length(data)),2)];
max_id = max(max(data_random));
A_random = sparse(data_random(:,1),data_random(:,2),1);
A_random(max_id,max_id) = 0 ; % Make sure the matrix is a square
A_random = A_random + A_random';
random_gcc = largestcomponent(A_random);
A_random_gcc = A_random(random_gcc,random_gcc); % Adjacency matrix of restricted to the GCC
dist_random_gcc = degree_distribution(A_random_gcc);
length(random_gcc)
%figure(2);
%h2 = plot(graph(A_random));
%layout(h2,'force','UseGravity',true)
%% configuration
dist = degree_distribution(A);
A_config = configuration(dist);
length(largestcomponent(A_config))
%figure(3);
%h3 = plot(graph(A_config));
%layout(h3,'force','UseGravity',true)
%% new method
k = [1:length(dist)]';
r = dist(:,2);
r = r/sum(r);
u = 0;
for i = 1:1000
    p = r./(1-u.^k);
    u = sum(k.*p.*u.^(k-1))/sum(k.*p);
end
s = 1/sum(p);

dist_full = [k,round(length(A)/s*p/sum(p))];
A_full = configuration(dist_full);
%figure(4);
%h4 = plot(graph(A_full));
%layout(h4,'force','UseGravity',true)

gcc = largestcomponent(A_full);
A_full_GCC = A_full(gcc,gcc);

dist_GCC = degree_distribution(A_full_GCC);
%%
figure(5); hold on;
fig1 = plot(dist(:,2));
fig2 = plot(dist_full(:,2),'^');
fig3 = plot(dist_GCC(:,2),'o');
fig4 = plot(dist_random_gcc(:,2),'o');

legend([fig1,fig2,fig3,fig4],'Original Network','Random Full','Random GCC','Shuffling GCC')

figure('Name', 'QQplot') ; hold on ;
qq1 = qqplot(dist(:, 2), [dist_full(:, 2), dist_GCC(:, 2), dist_random_gcc(:, 2)]) ;
legend('Random Full','Random GCC','Shuffling GCC');