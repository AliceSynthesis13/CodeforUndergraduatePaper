clear all;
clc;

linewidth = 3;
circlewidth = 30;
fontsize = 20;

Gamma = [0, 0];
a1 = [         0,    1];
a2 = [ sqrt(3)/2, -1/2];
a3 = [-sqrt(3)/2, -1/2];
b1 = a2 - a3;
b2 = a3 - a1;
b3 = a1 - a2;

figure(1);

quiver(0,0,a1(1),a1(2),0.96,'r','LineWidth',2); hold on;
quiver(0,0,a2(1),a2(2),0.96,'r','LineWidth',2); hold on;
quiver(0,0,a3(1),a3(2),0.96,'r','LineWidth',2); hold on;

quiver(a3(1),a3(2),b1(1),b1(2),0.96,'b','LineWidth',2); hold on;
quiver(a1(1),a1(2),b2(1),b2(2),0.96,'b','LineWidth',2); hold on;
quiver(a2(1),a2(2),b3(1),b3(2),0.96,'b','LineWidth',2); hold on;

scatter(0,0,'k','filled','SizeData',circlewidth); hold on;

text(a1(1)/3-0.2,a1(2)/3,'$a_1$','FontSize',24,'Interpreter','latex'); hold on;
text(a2(1)/3+0.12,a2(2)/3,'$a_2$','FontSize',24,'Interpreter','latex'); hold on;
text(a3(1)/3-0.24,a3(2)/3+0.05,'$a_3$','FontSize',24,'Interpreter','latex'); hold on;

text(-2*a1(1)/3-0.1,-2*a1(2)/3+0.06,'$b_1$','FontSize',24,'Interpreter','latex'); hold on;
text(-2*a2(1)/3-0.05,-2*a2(2)/3,'$b_2$','FontSize',24,'Interpreter','latex'); hold on;
text(-2*a3(1)/3-0.12,-2*a3(2)/3,'$b_3$','FontSize',24,'Interpreter','latex'); hold on;

axis equal;
axis off;
set(gcf,'unit','normalized','position',[0.3,0.3,0.3,0.3]);
saveas(gcf,'.\vector.jpg')
