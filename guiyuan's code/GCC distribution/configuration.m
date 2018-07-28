function A = configuration(dist)
number = dist(dist(:,2)>0,1:2);
x = [];
for i = 1:length(number)
    id = sum(number(1:i-1,2))+1:sum(number(1:i,2));
    temp = ones(number(i),1)*id;
    x = [x; temp(:)];
end
if mod(length(x),2)==1
    x = [1;x+1];
end
x = x(randperm(length(x)));
x = reshape(x,[length(x)/2,2]);
x(x(:,1)==x(:,2),:) = [];
A = sparse(x(:,1),x(:,2),1);
A = logical(A);
A = A + A';
A = logical(A);