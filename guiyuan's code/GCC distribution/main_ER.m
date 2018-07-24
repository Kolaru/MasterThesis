clear; clc;
max_k = 100;
k = 1:max_k;


c_list = [1:0.01:3]';
d_list = zeros(length(c_list),1);
u_list = zeros(length(c_list),1);
s_list = zeros(length(c_list),1);

for j = 1:length(c_list)
    
    c = c_list(j)
    
    r = c.^k*exp(-c)./factorial(k);
    
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
xlabel('c');
ylabel('average degree');

subplot(1,3,2); hold on;
plot(d_list,u_list);
xlabel('average degree');
ylabel('u');
grid on;

subplot(1,3,3); hold on;
plot(d_list,s_list);
xlabel('average degree');
ylabel('s');
grid on;
%%
c = 2;
r = c.^k*exp(-c)./factorial(k);
r = r/sum(r);
u = 0;
for i = 1:1000
    p = r./(1-u.^k);
    u = sum(k.*p.*u.^(k-1))/sum(k.*p);
end

figure(2);
pk = p/sum(p);
figure(2); hold on;
fig1 = plot(1:10,pk(1:10),'-o');
fig2 = plot(1:10,r(1:10),'-o');
xlabel('degree');
legend([fig1,fig2],'p(k)','r(k)')
legend boxoff;
set(gca,'FontSize',18)
box on;