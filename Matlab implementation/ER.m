function A = ER(N,c)
x = randi([1,N],round(1.01*c*N),2);
x( x(:,2)>=x(:,1),: )=[];
x(round(c*N/2)+1+3:end,:) = [];
A = sparse(x(:,1),x(:,2),1,N,N);
clear x;
A = A + A';
A = logical(A);