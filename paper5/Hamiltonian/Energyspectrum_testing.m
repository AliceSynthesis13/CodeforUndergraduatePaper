clear all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 以下画图程序只针对非开边界的能谱
% 这是看布里渊区内能带的

interval = 200; %interval 

kx1 = -pi/2 + zeros(1,interval);
ky1 = -pi/2 + 0:(pi/2/interval):(0-pi/2/interval);

kx2 = -pi/2:(pi/2/interval):(0-pi/2/interval);
ky2 = 0 + zeros(1,interval);

kx3 = 0:(pi/2/interval*sqrt(2)):pi/2;
ky3 = 0:(pi/2/interval*sqrt(2)):pi/2;

kx = [kx1, kx2, kx3];
ky = [ky1, ky2, ky3];

pointx = length(kx);
pointy = length(ky);

rows = (1:1:pointx)./pointx*(pi+pi/sqrt(2));

bands = 12;

energyspectrums = zeros(bands,pointx);

for x = 1 : pointx  
    k = [kx(x),ky(x)];
    E = Energyspectrum(k);
    energyspectrums(:,x) = E ;  
end

for i = 1 : bands
        plot(rows,energyspectrums(i,:),'k','linewidth',1) ; hold on
end


set(gca,'FontName','Times New Roman','FontSize',20);

% plot([0,0],[-4,4],'--b','linewidth',1);hold on
% plot([-4*pi/3,-4*pi/3],[-4,4],'--b','linewidth',1);hold on
% plot([-2*pi,2*pi/sqrt(3)],[0,0],'--r','linewidth',1);hold on
% xlim([-(2+sqrt(2))*pi,(2+sqrt(2))*pi]);

% set(gca,'XTickLabel',{});
xlim([0,(1+1/sqrt(2))*pi]);
xticks([0 pi/2 pi pi+pi/sqrt(2)]);
ylabel('\it{E}');

ax = gca;
ax.TickLabelInterpreter = 'latex';
xticklabels({'K','M','$\Gamma$','K'});

set(gcf,'unit','normalized','position',[0.1,0.1,0.4,0.4]);
saveas(gcf,'.\Energyspectrum_paper5.jpg')

% title(strcat(str,'-wave,',' SOC=',num2str(),'Fontsize',16));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear all;
% clc;
% 
% interval = 200;
% kx = -pi/2:pi/interval:pi/2;
% ky = -pi/2:pi/interval:pi/2;
% pointx = length(kx);
% pointy = length(ky);
% rows = 12;
% E1 = zeros(pointx,pointy,rows);
% E2 = zeros(pointx,pointy,rows);
% 
% for x = 1:pointx
%     for y = 1:pointy
%         E1(x,y,:) = Energyspectrum([ kx(x), ky(y)]);
%         E2(x,y,:) = Energyspectrum([ ky(y),-kx(x)]);
%     end
% end
% E = E1 - E2;
% 
% ftsize = 20;
% n = 6;
% figures = pcolor(kx,ky,E(:,:,n)); %,'ShowText','on'
% colorbar
% set(figures,'linestyle','none');
% maxv = max(max(E(:,:,n)));
% minv = min(min(E(:,:,n)));
% caxis([minv,maxv]);