clear all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 以下针对 x 开边界下的能谱
% y 方向周期性，可以用 ky 获得能谱
interval = 400; %interval

kx = 0:pi/interval:pi;

% type = 2 与 1 一样
type = 1;

n = 20;

str = strcat('.\Figure\OBC_M3_type',num2str(type),'.jpg');

if type == 1 || type == 2
    bands = 12*n;
elseif type == 3
    bands = 12*n-4;
elseif type == 4
    bands = 12*n-8;
end

rows = kx;

pointy = length(kx);    %Point Total

energyspectrums = zeros(bands,pointy);

for x = 1 : pointy
    k = [kx(x), 0];
    E = Energyspectrum(k,n,type);
    energyspectrums(:,x) = E ;
end

for i = 1 : bands
    plot(rows,energyspectrums(i,:),'k','linewidth',0.5) ; hold on
end

plot(rows,energyspectrums(bands/2-1,:),'r','linewidth',1) ; hold on
plot(rows,energyspectrums(bands/2  ,:),'r','linewidth',1) ; hold on
plot(rows,energyspectrums(bands/2+1,:),'r','linewidth',1) ; hold on
plot(rows,energyspectrums(bands/2+2,:),'r','linewidth',1) ; hold on




% plot([0,0],[-4,4],'--b','linewidth',1);hold on
plot([0,2*pi],[0,0],'--r','linewidth',1);hold on
xlim([0,pi]);
ylim([-2,2]);
xlabel('\it{k_x}');
ylabel('\it{E}');
xticks([0 pi/2 pi]);

ax = gca;
ax.TickLabelInterpreter = 'latex';
xticklabels({'0', '$\pi/2$','$\pi$'});

set(gca,'FontName','Times New Roman','FontSize',20);

% set(gca,'XTickLabel',{});

set(gcf,'unit','normalized','position',[0.2,0.2,0.4,0.5]);

saveas(gcf,str);