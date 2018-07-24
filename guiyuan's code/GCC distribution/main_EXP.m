clear; clc;
max_k = 100;
k = 1:max_k;


c_list = [1:0.001:3]';
d_list = zeros(length(c_list),1);
u_list = zeros(length(c_list),1);
s_list = zeros(length(c_list),1);

for j = 1:length(c_list)
    
    c = c_list(j)
    
    r = (1-1./c).^(k-1)./c;
    
    r = r/sum(r);
    
    d_list(j) = sum(k.*r);
    
    
    u = 0;
    for i = 1:1000
        p = r./(1-u.^k);
        u = sum(k.*p.*u.^(k-1))/sum(k.*p);
    end
    
    u_list(j) = u;
    s_list(j) = 1/sum(p);
    
end

figure(1); hold on;
subplot(1,3,1); hold on;
plot(c_list,d_list);
grid on;

subplot(1,3,2); hold on;
plot(d_list,u_list);
grid on;

subplot(1,3,3); hold on;
plot(d_list,s_list);
grid on;