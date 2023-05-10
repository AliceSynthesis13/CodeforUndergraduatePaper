clear all;
clc;

global Nx;
global Ny;
global Nmovs;
global Nmov;

% 这里的 Nx 是函数里的 Ny，也就是竖直方向的原子数目
% 这里的 Ny 是函数里的 Nx，也就是水平方向的原子数目
Nx = 21;
Ny = 20;

Nmovs = 1; % Nmovs为起始添加sublattice的位置
Nmov  = 0; % Nmov为添加的sublattice的个数, Nmov < Ny以便观察到现象

% type = 1 对应的是边界均为 A-B 型，哈密顿量不需要改变
% type = 2 对应的是边界均为 A-A 型，哈密顿量需要去掉最下面的 A 层和最上面的 B 层

type = 2;

N0 = type - 1;
%首先导出哈密顿量并对角化
% H = Real_Hamiltonian(Nx,Ny,type);
% [V,D] = eig(kron(H,eye(2)));                         % V 是本征向量，D 是本征值
% [Y,I] = sort(diag(D),'descend');        % 降序，Y 是已经排序好的数组，而 I 是 Y 每个数在原来数组（矩阵）中的位置
% min(abs(Y));

r=3.5;
pi=2*3.14159;
l=3;
csize=100;
sensitive=1;
linewidth=0.5;
circlewidth=0.3;

spectrum = 1; % 是否输出能谱(0为不输出,1为输出)
plane = 0;    % 是否输出平面图
lattice = 0;  % 是否输出晶格图(在plane=1的情况下)

matrix = Real_Hamiltonian(Ny,Nx,type);
matrix = kron(matrix,[1,0;0,-1]);
[vet,value]=eig(matrix);

[Y,I] = sort(diag(value),'descend');        % 降序，Y 是已经排序好的数组，而 I 是 Y 每个数在原来数组（矩阵）中的位置
L = length(Y);
min(abs(Y));
Value = Y(L/2-5:L/2+6)'

if spectrum==1
    figure(2);
    scatter(-5.5:2:25.5,-Y(length(Y)/2-15:2:length(Y)/2+16),'linewidth',1);hold on;
    scatter(8.5:2:11.5,-Y(length(Y)/2-1:2:length(Y)/2+2),'red','linewidth',2);
    hold on;
    set(gca,'FontName','Times New Roman','FontSize',20);
    xticklabels({''})
    yticklabels({''})
    ylim([-4e-1,4e-1]);
    yticks([-4e-1 0 4e-1]);
    box on;

    ay = gca;
    ay.TickLabelInterpreter = 'latex';
    yticklabels({'-0.4', '0','0.4'});
    ylabel('\it{E}');
    set(gcf,'unit','normalized','position',[0.4,0.4,0.3,0.36]);
    set(gca, 'LooseInset', [0,0,0,0]);
    set(0,'defaultfigurecolor','w');
    % saveas(gcf,'.\Figure\R_OBC_M3_DoS_type3.jpg');
    export_fig('.\Figure\R_OBC_M20_E_type3.jpg');
end

Phi = [];
minex = min(abs(Y))
for i = 1:length(value(1,:))
    if abs(value(i,i)) < 1.0001*minex
        Phi = [Phi; vet(:,i)'];
    end
end

U = zeros(1,length(Phi(1,:)));
for i = 1:length(Phi(1,:))
    U(i) = Phi(:,i)'*Phi(:,i);
end

S = zeros(1,length(Phi(1,:))/4);
for i = 1:length(Phi(1,:))/4
    S(i) = sum(U((4*(i-1)+1):(4*(i-1)+4)));
end

% S = S./sum(S);

tot = 0;
if plane == 1
    y0 = 0;
    Max = 0;
    figure(1);
    for i=1:1:Nx
        if (mod(N0+1,4)==mod(i,4)||mod(N0+3,4)==mod(i,4))
            y0=y0-l;
        end
        if (mod(N0+2,4)==mod(i,4)||mod(N0+4,4)==mod(i,4))
            y0=y0-l/2;
        end
        for j=1:1:Ny
            if (mod(N0+1,4)==mod(i,4)||mod(N0+4,4)==mod(i,4))
                x(Ny*(i-1)+j,1)=l*sqrt(3)*(j-1/2);
                y(Ny*(i-1)+j,1)=y0;
            end
            if (mod(N0+2,4)==mod(i,4)||mod(N0+3,4)==mod(i,4))
                x(Ny*(i-1)+j,1)=l*sqrt(3)*(j-1);
                y(Ny*(i-1)+j,1)=y0;
            end
            bohanshu=0;
            for k=1:1:4
                %               for k2=1:1:2%加磁场
                for k2=0:1:3%正常情况
                    bohanshu=bohanshu+abs(vet(posi(i,j)+k,2*(Nx*Ny+Nmov)-1+k2))^2;
                end
            end
            tot = tot + bohanshu;
            C(Ny*(i-1)+j,2)=sensitive*bohanshu;
            C(Ny*(i-1)+j,3)=C(Ny*(i-1)+j,2);
            if  C(Ny*(i-1)+j,2)>Max
                Max=C(Ny*(i-1)+j,2);
            end
            % 画晶格形状:（sublattice之间的连线）
            if mod(N0+1,4)==mod(i,4)
                if i>1
                    line([x(Ny*(i-1)+j),x(Ny*(i-1)+j)],[y0,y0+l],'Color','black','LineWidth',linewidth);
                    hold on;
                end
            elseif mod(N0+2,4)==mod(i,4)
                if i>1
                    line([x(Ny*(i-1)+j),x(Ny*(i-1)+j)+l/2*sqrt(3)],[y0,y0+l/2],'Color','black','LineWidth',linewidth);
                    hold on;
                    if j>1
                        line([x(Ny*(i-1)+j),x(Ny*(i-1)+j)-l/2*sqrt(3)],[y0,y0+l/2],'Color','black','LineWidth',linewidth);
                        hold on;
                    end
                end
            elseif mod(N0+3,4)==mod(i,4)
                if i>1
                    line([x(Ny*(i-1)+j),x(Ny*(i-1)+j)],[y0,y0+l],'Color','black','LineWidth',linewidth);
                    hold on;
                end
            elseif mod(N0+4,4)==mod(i,4)
                if i>1
                    line([x(Ny*(i-1)+j),x(Ny*(i-1)+j)-l/2*sqrt(3)],[y0,y0+l/2],'Color','black','LineWidth',linewidth);
                    hold on;
                    if j<Ny
                        line([x(Ny*(i-1)+j),x(Ny*(i-1)+j)+l/2*sqrt(3)],[y0,y0+l/2],'Color','black','LineWidth',linewidth);
                        hold on;
                    end
                end
            end
        end
    end
    y0=0;
    for j=1:1:Nmov % 对最上面一行的操作
        if N0==0
            x(Nx*Ny+j,1)=l*sqrt(3)*(j-3/2+Nmovs);
            y(Nx*Ny+j,1)=y0;
            % 画晶格形状:
            line([x(Nx*Ny+j),x(Nx*Ny+j)],[y0,y0-l],'Color','black','LineWidth',linewidth);
            hold on;
        end
        if N0==1
            x(Nx*Ny+j,1)=l*sqrt(3)*(j-2+Nmovs);
            y(Nx*Ny+j,1)=y0;
            % 画晶格形状:
            line([x(Nx*Ny+j),x(Nx*Ny+j)-l/2*sqrt(3)],[y0,y0-l/2],'Color','black','LineWidth',linewidth);
            hold on;
            if j<Ny
                line([x(Nx*Ny+j),x(Nx*Ny+j)+l/2*sqrt(3)],[y0,y0-l/2],'Color','black','LineWidth',linewidth);
                hold on;
            end
        end
        bohanshu=0;
        for k=1:1:4
            %       for k2=1:1:2%加磁场
            for k2=0:1:3%正常情况
                bohanshu=bohanshu+abs(vet(4*(Nx*Ny+j-1)+k,2*(Nx*Ny+Nmov)-1+k2))^2;
            end
        end
        tot = tot + bohanshu;
        C(Nx*Ny+j,2)=sensitive*bohanshu;
        C(Nx*Ny+j,3)=C(Nx*Ny+j,2);
        if  C(Nx*Ny+j,2)>Max
            Max=C(Nx*Ny+j,2);
        end
    end

    C(:,2) = S;
    C(:,3) = S;
    Max = max(S);
    C(:,2)=1/Max*(Max-C(:,2));
    C(:,3)=1/Max*(Max-C(:,3));
    C(:,1)=1; % 红色分布
    % C(:,1)=C(:,3);
    % C(:,2)=1-(255-26)/255*(C(:,2)/Max);
    % C(:,3)=1-(255-139)/255*(C(:,3)/Max);
    % C(:,1)=1-(255-85)/255*(C(:,1)/Max); % 紫色分布
    % C(:,2)=1;
    % C(:,1)=1/Max*(Max-C(:,3));
    % C(:,3)=1; % Cyan色分布
    scatter(x,y,csize,'k','LineWidth',circlewidth);
    hold on;
    if lattice~=1
        scatter(x,y,csize,C,'filled');
    else
        for i=1:1:Nx
            if  (mod(N0+1,4)==mod(i,4)||mod(N0+3,4)==mod(i,4))
                scatter(x((i-1)*Ny+1:1:i*Ny),y((i-1)*Ny+1:1:i*Ny),csize,'r','filled');
                hold on;
            elseif (mod(N0+2,4)==mod(i,4)||mod(N0+4,4)==mod(i,4))
                scatter(x((i-1)*Ny+1:1:i*Ny),y((i-1)*Ny+1:1:i*Ny),csize,'b','filled');
                hold on;
            end
        end
        scatter(x(Nx*Ny+1:1:Nx*Ny+Nmov),y(Nx*Ny+1:1:Nx*Ny+Nmov),csize,'b','filled');
        hold on;
    end
    axis equal;
    axis off;
    set(gcf,'unit','normalized','position',[0.1,0.1,0.8,0.7]);
    % set(gcf,'unit','normalized','position',[0.1,0.1,0.4,0.35]);
    set(gca, 'LooseInset', [0,0,0,0]);
    set(0,'defaultfigurecolor','w');
    % saveas(gcf,'.\Figure\R_OBC_M3_DoS_type3.jpg');
    export_fig('.\Figure\R_OBC_M20_DoS_type3.jpg');
    % export_fig('.\Figure\Lattice3.jpg');
end

%% 自定义函数

% 给定一个坐标 (i,j)，返回矩阵中的位置，针对矩形部分
% Nx 是从上到下的层数，Ny 是从左到右的行数
% i 是上下方向的标记，j 是左右方向的标记
function P=posi(i,j)
    global Nx; global Ny;
    P=4*(Ny*(i-1)+(j-1));
end

% 给定一个坐标 (i,j)，返回矩阵中的位置，针对额外添加在最上面的那一层部分
function P1=posim(j)
    global Ny; global Nx;global Nmovs;
    if j>=Nmovs
        P1=4*Nx*Ny+4*(j-Nmovs);
    end
    if j<Nmovs
        P1=4*Nx*Ny+4*(j+Ny-Nmovs);
    end
end

function M=mod1(i)
    global Ny;
    M=i;
    if i>Ny
        M=i-Ny;
    end
    if i<1
        M=i+Ny;
    end
end