clear all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 以下针对 x 开边界下的能谱
% y 方向周期性，可以用 ky 获得能谱
% 这是计算 Fig2 所采用的程序
interval = 400; %interval
n = 40;

% % type = 1 代表 armchair edges
% % type = 2 代表 zigzag edges（连续的锯齿）（双闭）KM 模型经典边界
% % type = 3 代表像栏杆一样的（双开）
% % type = 4 （单开单闭）

type = 4;

if type == 1 || type == 2
    bands = 8*n;
elseif type == 3
    bands = 8*n-8;
elseif type == 4
    bands = 8*n-4;
end

str = strcat('.\Figure\OBC_M10_type',num2str(type),'.jpg');

if type == 2 || type == 3 || type == 4
    kx = 0:pi/interval/sqrt(3):2*pi/sqrt(3);
    rows = kx;
    pointy = length(kx);    % Point Total

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


    set(gca,'FontSize',20,'FontName','Arial');

    % plot([0,0],[-4,4],'--b','linewidth',1);hold on
    plot([0,2*pi/sqrt(3)],[0,0],'--r','linewidth',1);hold on
    xlim([0,2*pi/sqrt(3)]);
    ylim([-2,2]);
    xlabel('\it{k_x}');
    ylabel('\it{E}');
    xticks([0 pi/sqrt(3) 2*pi/sqrt(3)]);

    ax = gca;
    ax.TickLabelInterpreter = 'latex';
    xticklabels({'0', '$\pi/\sqrt{3}$','$2\pi/\sqrt{3}$'});

    set(gca,'FontName','Times New Roman','FontSize',20);
    % set(gca,'XTickLabel',{});
    set(gcf,'unit','normalized','position',[0.2,0.2,0.6,0.6]);
    saveas(gcf,str);

elseif type == 1
    kx = -pi/3:pi/interval/3:pi/3;
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


    set(gca,'FontSize',20,'FontName','Arial');

    % plot([0,0],[-4,4],'--b','linewidth',1);hold on
    plot([-pi/3,pi/3],[0,0],'--r','linewidth',1);hold on
    xlim([-pi/3,pi/3]);
    ylim([-2,2]);
    xlabel('\it{k_x}');
    ylabel('\it{E}');
    xticks([-pi/3 0 pi/3]);

    ax = gca;
    ax.TickLabelInterpreter = 'latex';
    xticklabels({'$-\pi/3$', '0','$\pi/3$'});

    set(gca,'FontName','Times New Roman','FontSize',20);
    % set(gca,'XTickLabel',{});
    set(gcf,'unit','normalized','position',[0.2,0.2,0.6,0.6]);
    
    saveas(gcf,str);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
