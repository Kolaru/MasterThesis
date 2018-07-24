clear; clc;
A = ER(1000,1.5);
b = largestcomponent(A);
B = A(b,b);

dist_A = degree_distribution(A);
dist_B = degree_distribution(B);

figure(1); hold on;
set(gcf, 'position', [100 100 1700 700]);
subplot(1,2,1); hold on;
h = plot(graph(A));
layout(h,'force','UseGravity',true)
box on;
set(gca,'LineWidth',1.5);
%title('a','FontSize',20,'FontWeight','bold');
set(gca,'xtick',[])
set(gca,'ytick',[])
axis off;
%%
ylim=get(gca,'ylim');
xlim=get(gca,'xlim');
delta = (ylim(2) - ylim(1))*0.05;
text(xlim(1),ylim(2)+delta,'a','FontSize',20,'FontWeight','bold')


subplot(1,2,2); hold on;
fig1 = plot(dist_A(:,1),dist_A(:,3)/100,'-o','LineWidth',2);
fig2 = plot(dist_A(2:end,1),dist_A(2:end,3)/sum(dist_A(2:end,3)),'-o','LineWidth',2);
fig3 = plot(dist_B(:,1),dist_B(:,3)/100,'-o','LineWidth',2);
legend([fig1,fig2,fig3],'Original Distribution','Original Distribution(d > 0)','GCC Distributions')
legend boxoff;
box on;
set(gca,'LineWidth',1.5);
set(gca,'FontSize',18);
xlabel('degree')
ylabel('P')
%title('b','FontSize',20,'FontWeight','bold');

ylim=get(gca,'ylim');
xlim=get(gca,'xlim');
delta = (ylim(2) - ylim(1))*0.05;
text(xlim(1),ylim(2)+delta,'b','FontSize',20,'FontWeight','bold')
