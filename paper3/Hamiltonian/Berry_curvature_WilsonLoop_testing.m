% function[Chern] = Berry_curvature_WilsonLoop_testing(interval)
%%% 本程序利用 Wilson Loop 的方式计算逐点的 Berry curvature
%%% 本程序采用与之前程序一脉相承的思路 用于计算有简并能谱的 Berry curvature

clear all; clc

interval = 2/200;
dkx = interval*4/3*pi;       %interval
dky = interval*2/sqrt(3)*pi; %interval

kix = -4*pi/3;   %initial range
kfx =  4*pi/3;    %final range
kiy = -2*pi/sqrt(3);
kfy =  2*pi/sqrt(3);

kx = kix:dkx:kfx;
ky = kiy:dky:kfy;

pointx = length(kx);    %Point Total
pointy = length(ky);    %Point Total


rows = 4;
Berry = zeros(max(pointx,pointy),max(pointx,pointy),rows);

for y = 1:pointy
    
    if y>=1 && y<pointy/2
        y3=-(y*dky+kiy)/(dkx*sqrt(3)) ;
        y4= (y*dky+kiy)/(dkx*sqrt(3)) - 2*kix/(dkx);
        y3=int64(y3);
        y4=int64(y4);
        for x = y3+1:y4+1
            R=[ky(y),kx(x)];
            Berry(y,x,:) = Berry_curvature_WilsonLoop2(R,dkx,dky,rows);
        end
    end
    
    if y>pointy/2 && y<pointy
        y1= (y*dky+kiy)/(dkx*sqrt(3)) ;
        y2=-(y*dky+kiy)/(dkx*sqrt(3)) - 2*kix/(dkx);
        y1=int64(y1);
        y2=int64(y2);
        for x = y1+1:y2+1
            R=[ky(y),kx(x)];
            Berry(y,x,:) = Berry_curvature_WilsonLoop2(R,dkx,dky,rows);
        end
    end
end

Berry = flip(Berry,3);

Chern = zeros(1,rows);
for m = 1 : rows
    for y = 1:pointy
        for x = 1:pointx
            Chern(m) = Chern(m) + 1/(2*pi) * Berry(y,x,m) * dkx * dky; % SOC only时要乘4！！！
        end
    end
end


%画图
ftsize = 20;

figures = pcolor(ky,kx,Berry(:,:,2)); %,'ShowText','on'
colorbar
% xlabel('$k_x$','Fontsize',ftsize,'Interpreter','latex');
% ylabel('$k_y$','Fontsize',ftsize,'Interpreter','latex');
% set(gca,'linewidth',2,'fontsize',ftsize,'fontname','Times');
% 
% set(gca,'XTick',{});
% set(gca,'YTick',{});
% 
% text(-2*pi/3-0.3,-pi/sqrt(3)-0.2,'$-\frac{2}{3}\pi$','fontsize',30,'Interpreter','latex');
% text(-1*pi/3-0.3,-pi/sqrt(3)-0.2,'$-\frac{1}{3}\pi$','fontsize',30,'Interpreter','latex');
% %text(-0.1       ,-pi/sqrt(3)-0.2,'0','fontsize',30);
% text(+1*pi/3-0.2,-pi/sqrt(3)-0.2,'$\frac{1}{3}\pi$','fontsize',30,'Interpreter','latex');
% text(+2*pi/3-0.2,-pi/sqrt(3)-0.2,'$\frac{2}{3}\pi$','fontsize',30,'Interpreter','latex');
% 
% text(-2*pi/3-0.7,-1*pi/sqrt(3),'$-\frac{1}{\sqrt{3}}\pi$','fontsize',30,'Interpreter','latex');
% text(-2*pi/3-0.7,-0.5*pi/sqrt(3),'$-\frac{1}{2\sqrt{3}}\pi$','fontsize',30,'Interpreter','latex');
% %text(-2*pi/3-0.8,-0*pi/sqrt(3),'0','fontsize',30,'Interpreter','latex');
% text(-2*pi/3-0.6,+0.5*pi/sqrt(3),'$\frac{1}{2\sqrt{3}}\pi$','fontsize',30,'Interpreter','latex');
% text(-2*pi/3-0.5,+1*pi/sqrt(3),'$\frac{1}{\sqrt{3}}\pi$','fontsize',30,'Interpreter','latex');
% 
% % xlim([-2*pi/3,2*pi/3]);
% % ylim([-1*pi/sqrt(3),1*pi/sqrt(3)]);
% xlim([-2*pi/3,2*pi/3]);
% ylim([-1*pi/sqrt(3),1*pi/sqrt(3)]);% for SOC only
caxis([-1,1]);
set(figures,'linestyle','none');
% set(gcf,'unit','normalized','position',[0.2,0.2,0.52,0.64]);








