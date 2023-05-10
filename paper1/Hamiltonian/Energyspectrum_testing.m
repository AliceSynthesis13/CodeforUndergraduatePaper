clear all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 以下画图程序只针对非开边界的能谱
% 这是看布里渊区内能带的

interval = 1000; %interval 

kx1 = 0:pi/interval:2*sqrt(3)*pi/3;
ky1 = zeros(1,length(kx1));

ky2 = 0:pi/interval:2*pi/3;
kx2 = zeros(1,length(ky2)) + 2*sqrt(3)*pi/3;

kx = [kx1, kx2];
ky = [ky1, ky2];

rows = (1:1:length(kx))./length(kx).*(2*sqrt(3)*pi/3+2*pi/3);

pointx = length(kx);    %Point Total
pointy = length(ky);    %Point Total

bands = 4;

energyspectrums = zeros(bands,pointx);

for x = 1 : length(kx)   
    k = [kx(x),ky(x)];  
    E = Energyspectrum(k);
    energyspectrums(:,x) = E;    
end

for i = 1 : bands
        plot(rows,energyspectrums(i,:),'k','linewidth',2) ; hold on;
end


set(gca,'FontSize',20,'FontName','Arial');


plot([4*sqrt(3)*pi/9,4*sqrt(3)*pi/9],[-4,4],'--b','linewidth',1);hold on;
plot([2*sqrt(3)*pi/3,2*sqrt(3)*pi/3],[-4,4],'--b','linewidth',1);hold on;
% Fermi surface
plot([0,(2*sqrt(3)*pi/3+2*pi/3)],[0,0],'--r','linewidth',1);hold on;
xlim([0,(2*sqrt(3)*pi/3+2*pi/3)]);

xticks([0 4*sqrt(3)*pi/9 2*sqrt(3)*pi/3 (2*sqrt(3)*pi/3+2*pi/3)]);

ax = gca;
ax.TickLabelInterpreter = 'latex';
xticklabels({'$\Gamma$', 'K','M','$\Gamma$'});

set(gcf,'unit','normalized','position',[0.1,0.1,0.4,0.4]);
% saveas(gcf,'.\Energyspectrum_paper2.jpg')
% title(strcat(str,'-wave,',' SOC=',num2str(),'Fontsize',16));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 这是看 y 方向开边界的
% interval = 400; %interval 
% 
% kx = -pi:pi/interval:pi;
% 
% rows = kx;
% 
% pointy = length(kx);    %Point Total
% 
% bands = 4;
% 
% energyspectrums = zeros(bands,pointy);
% 
% for x = 1 : pointy
%     
%     k = [kx(x), 0];
%    
%     E = Energyspectrum(k);
%     
%     energyspectrums(:,x) = E ;
%     
% end
% 
% 
% for i = 1 : bands
%  
%         plot(rows,energyspectrums(i,:),'k','linewidth',2) ; hold on
% end
% 
% 
% set(gca,'FontSize',20,'FontName','Arial');
% 
% 
% plot([0,0],[-4,4],'--b','linewidth',1);hold on
% plot([-pi,pi],[0,0],'--r','linewidth',1);hold on
% xlim([-pi,pi]);
% ylim([-4,4]);
% 
% % set(gca,'XTickLabel',{});
% 
% set(gcf,'unit','normalized','position',[0.2,0.2,0.52,0.74]);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % interval = 2/600; %interval
% % 
% % m = pi*sqrt(3);
% % 
% % dkx = interval*m;       %interval
% % dky = interval*m*2/sqrt(3); %interval
% % 
% % kix = -m;    %initial range
% % kfx =  m;    %final range
% % kiy = -m*2/sqrt(3);
% % kfy =  m*2/sqrt(3);
% % 
% % kx = kix:dkx:kfx;
% % ky = kiy:dky:kfy;
% % 
% % pointx = length(kx);    %Point Total
% % pointy = length(ky);    %Point Tota
% % 
% % E = zeros(pointx,pointy,4);
% % 
% % for x = 1:pointx
% %     for y = 1:pointy
% %         E(y,x,:) = Energyspectrum([kx(x),ky(y)]);
% %     end
% % end
% % 
% % ftsize = 20;
% % 
% % figure = pcolor(kx,ky,E(:,:,2)) %,'ShowText','on'
% % colorbar;
% % xlabel('k_x','Fontsize',ftsize);
% % ylabel('k_y','Fontsize',ftsize);
% % set(gca,'linewidth',2,'fontsize',ftsize,'fontname','Times');
% % % set(gca,'XTick',-1:0.5:1);
% % % set(gca,'YTick',-1:0.5:1);
% % 
% % caxis([0,3]);
% % set(figure,'linestyle','none');
% % 
% % 
% % 
