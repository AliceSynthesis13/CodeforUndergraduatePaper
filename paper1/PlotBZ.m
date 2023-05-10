% 这个文件是用来画布里渊区的示意图
clear all;
clc;

linewidth = 3;
circlewidth = 100;
fontsize = 20;

K1 = [-2*pi/3/sqrt(3), -2*pi/3];
K2 = [ 2*pi/3/sqrt(3), -2*pi/3];
K3 = [ 4*pi/3/sqrt(3),       0];
K4 = [ 2*pi/3/sqrt(3),  2*pi/3];
K5 = [-2*pi/3/sqrt(3),  2*pi/3];
K6 = [-4*pi/3/sqrt(3),       0];

M1 = [0, -2*pi/3];

figure(1);

line([K1(1),K2(1)],[K1(2),K2(2)],'Color','black','LineWidth',linewidth); hold on;
line([K2(1),K3(1)],[K2(2),K3(2)],'Color','black','LineWidth',linewidth); hold on;
line([K3(1),K4(1)],[K3(2),K4(2)],'Color','black','LineWidth',linewidth); hold on;
line([K4(1),K5(1)],[K4(2),K5(2)],'Color','black','LineWidth',linewidth); hold on;
line([K5(1),K6(1)],[K5(2),K6(2)],'Color','black','LineWidth',linewidth); hold on;
line([K6(1),K1(1)],[K6(2),K1(2)],'Color','black','LineWidth',linewidth); hold on;

scatter(K1(1),K1(2),'r','filled','SizeData',circlewidth); hold on;
text(K1(1)-0.4,K1(2)-0.4,'K','FontSize',30,'Interpreter','latex'); hold on;

scatter(0,0,'r','filled','SizeData',circlewidth); hold on;
text(-0.5,0,'$\Gamma$','FontSize',30,'Interpreter','latex'); hold on;

scatter(0,-2*pi/3,'r','filled','SizeData',circlewidth); hold on;
text(M1(1)-0.2,M1(2)-0.4,'M','FontSize',30,'Interpreter','latex'); hold on;

quiver(0,0,K1(1),K1(2),1,'m','LineWidth',2); hold on;
quiver(K1(1),K1(2),2*pi/3/sqrt(3),0,1,'m','LineWidth',2); hold on;
quiver(M1(1),M1(2),0,2*pi/3,1,'m','LineWidth',2);

scatter(K3(1),K3(2),'r','filled','SizeData',circlewidth); hold on;
text(K3(1)-0.10,K3(2)-0.4,'$\frac{4\sqrt{3}\pi}{9}$','FontSize',30,'Interpreter','latex'); hold on;

axis equal;
axis off;
% set(gcf,'unit','normalized','position',[0.3,0.3,0.4,0.4]);
saveas(gcf,'.\BZ1.jpg')
