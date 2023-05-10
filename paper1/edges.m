clear all;
clc;

linewidth = 3;

Gamma1 = [  0,        0];
Gamma2 = [3/2,sqrt(3)/2];
Gamma3 = [  3,        0];
Gamma4 = [  3,  sqrt(3)];

K1 = [   1,          0];
K2 = [ 1/2,  sqrt(3)/2];
K3 = [-1/2,  sqrt(3)/2];
K4 = [  -1,          0];
K5 = [-1/2, -sqrt(3)/2];
K6 = [ 1/2, -sqrt(3)/2];

Point1 = [K1;K2;K3;K4;K5;K6];
Point2 = Point1 + kron(ones(6,1),Gamma2);
Point3 = Point1 + kron(ones(6,1),Gamma3);
Point4 = Point1 + kron(ones(6,1),Gamma4);

% Black Line
line([Point1(1,1),Point1(2,1)],[Point1(1,2),Point1(2,2)],'Color','black','LineWidth',linewidth); hold on;
line([Point1(2,1),Point1(3,1)],[Point1(2,2),Point1(3,2)],'Color','black','LineWidth',linewidth); hold on;
line([Point1(3,1),Point1(4,1)],[Point1(3,2),Point1(4,2)],'Color','black','LineWidth',linewidth); hold on;
line([Point1(4,1),Point1(5,1)],[Point1(4,2),Point1(5,2)],'Color','black','LineWidth',linewidth); hold on;
line([Point1(5,1),Point1(6,1)],[Point1(5,2),Point1(6,2)],'Color','black','LineWidth',linewidth); hold on;
line([Point1(6,1),Point1(1,1)],[Point1(6,2),Point1(1,2)],'Color','black','LineWidth',linewidth); hold on;

line([Point2(1,1),Point2(2,1)],[Point2(1,2),Point2(2,2)],'Color','black','LineWidth',linewidth); hold on;
line([Point2(2,1),Point2(3,1)],[Point2(2,2),Point2(3,2)],'Color','black','LineWidth',linewidth); hold on;
line([Point2(3,1),Point2(4,1)],[Point2(3,2),Point2(4,2)],'Color','black','LineWidth',linewidth); hold on;
line([Point2(4,1),Point2(5,1)],[Point2(4,2),Point2(5,2)],'Color','black','LineWidth',linewidth); hold on;
line([Point2(5,1),Point2(6,1)],[Point2(5,2),Point2(6,2)],'Color','black','LineWidth',linewidth); hold on;
line([Point2(6,1),Point2(1,1)],[Point2(6,2),Point2(1,2)],'Color','black','LineWidth',linewidth); hold on;

line([Point3(1,1),Point3(2,1)],[Point3(1,2),Point3(2,2)],'Color','black','LineWidth',linewidth); hold on;
line([Point3(2,1),Point3(3,1)],[Point3(2,2),Point3(3,2)],'Color','black','LineWidth',linewidth); hold on;
line([Point3(3,1),Point3(4,1)],[Point3(3,2),Point3(4,2)],'Color','black','LineWidth',linewidth); hold on;
line([Point3(4,1),Point3(5,1)],[Point3(4,2),Point3(5,2)],'Color','black','LineWidth',linewidth); hold on;
line([Point3(5,1),Point3(6,1)],[Point3(5,2),Point3(6,2)],'Color','black','LineWidth',linewidth); hold on;
line([Point3(6,1),Point3(1,1)],[Point3(6,2),Point3(1,2)],'Color','black','LineWidth',linewidth); hold on;

line([Point4(1,1),Point4(2,1)],[Point4(1,2),Point4(2,2)],'Color','black','LineWidth',linewidth); hold on;
line([Point4(2,1),Point4(3,1)],[Point4(2,2),Point4(3,2)],'Color','black','LineWidth',linewidth); hold on;
line([Point4(3,1),Point4(4,1)],[Point4(3,2),Point4(4,2)],'Color','black','LineWidth',linewidth); hold on;
line([Point4(4,1),Point4(5,1)],[Point4(4,2),Point4(5,2)],'Color','black','LineWidth',linewidth); hold on;
line([Point4(5,1),Point4(6,1)],[Point4(5,2),Point4(6,2)],'Color','black','LineWidth',linewidth); hold on;
line([Point4(6,1),Point4(1,1)],[Point4(6,2),Point4(1,2)],'Color','black','LineWidth',linewidth); hold on;

% Red Line
line([Point1(4,1),Point1(5,1)],[Point1(4,2),Point1(5,2)],'Color','red','LineWidth',linewidth); hold on;
line([Point1(5,1),Point1(6,1)],[Point1(5,2),Point1(6,2)],'Color','red','LineWidth',linewidth); hold on;
line([Point1(6,1),Point1(1,1)],[Point1(6,2),Point1(1,2)],'Color','red','LineWidth',linewidth); hold on;

line([Point2(5,1),Point2(6,1)],[Point2(5,2),Point2(6,2)],'Color','red','LineWidth',linewidth); hold on;

line([Point3(4,1),Point3(5,1)],[Point3(4,2),Point3(5,2)],'Color','red','LineWidth',linewidth); hold on;
line([Point3(5,1),Point3(6,1)],[Point3(5,2),Point3(6,2)],'Color','red','LineWidth',linewidth); hold on;
line([Point3(6,1),Point3(1,1)],[Point3(6,2),Point3(1,2)],'Color','red','LineWidth',linewidth); hold on;

% Blue Line
line([Point3(1,1),Point3(2,1)],[Point3(1,2),Point3(2,2)],'Color','blue','LineWidth',linewidth); hold on;
% line([Point3(6,1),Point3(1,1)],[Point3(6,2),Point3(1,2)],'Color','blue','LineWidth',linewidth); hold on;

line([Point4(1,1),Point4(2,1)],[Point4(1,2),Point4(2,2)],'Color','blue','LineWidth',linewidth); hold on;
line([Point4(6,1),Point4(1,1)],[Point4(6,2),Point4(1,2)],'Color','blue','LineWidth',linewidth); hold on;

% 图注
text(1/2,-1.0,'armchair edges','FontSize',30,'Interpreter','latex'); hold on;
text(2.6,2.8,'zigzag edges','FontSize',30,'Interpreter','latex'); hold on;

axis equal;
axis off;
saveas(gcf,'.\Edges.jpg')
