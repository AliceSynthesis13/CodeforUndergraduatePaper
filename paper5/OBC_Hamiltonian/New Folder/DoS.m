clear all; clc;

t = 1;
lambda_SO = 0.4;
mu = 0.5;
Delta = 0.5;

Nx = 30;
Ny = 10;
N_begin = 7;
N_add = 20;
freedom = 4;

% 这里分两种画图，一种是示意图，一种是态密度图
% lattice = 1 代表示意图；lattice = 0 代表态相关的图；
% spectrum = 1 代表输出能量分布图 E ，spectrum = 0 代表输出态密度分布图 DoS
lattice = 0;
spectrum = 0;

% type 的输入值有 1、2、3、4 四种，对应 AB-C \ C-AB \ AB-AB \ C-C 四种类型
type = 4;

str1 = strcat('.\Figure\R_OBC_M8_E_type',num2str(type),'.jpg');
str2 = strcat('.\Figure\R_OBC_M8_DoS_type',num2str(type),'.jpg');
% str2 = '.\Figure\type4.jpg';

%首先导出哈密顿量并对角化
H = Real_Hamiltonian(Nx,Ny,type,t,lambda_SO,mu,Delta,N_begin,N_add);
[V,D] = eig(H);                         % V 是本征向量，D 是本征值
[Y,I] = sort(diag(D),'descend');        % 降序，Y 是已经排序好的数组，而 I 是 Y 每个数在原来数组（矩阵）中的位置
L = length(Y);
min(abs(Y));
Value = Y(L/2-5:L/2+6)'

if spectrum == 1
    figure(2);
    scatter(-9.5:1:29.5,-Y(length(Y)/2-19:length(Y)/2+20),'linewidth',1);hold on;
    % scatter(4.5:1:15.5,-Y(length(Y)/2-5:length(Y)/2+6),'red','linewidth',2);
    scatter(8.5:1:11.5,-Y(length(Y)/2-1:length(Y)/2+2),'red','linewidth',2);
    hold on;
    set(gca,'FontName','Times New Roman','FontSize',20);
    xticklabels({''})
    yticklabels({''})
    rangev = 1.5e-1;
    ylim([-rangev,rangev]);
    yticks([-rangev 0 rangev]);
    box on;

    ay = gca;
    ay.TickLabelInterpreter = 'latex';
    nv = strcat('-',num2str(rangev));
    pv = strcat(num2str(rangev));
    yticklabels({nv, '0',pv});
    ylabel('\it{E}');
    set(gcf,'unit','normalized','position',[0.4,0.4,0.3,0.36]);
    saveas(gcf,str1);

elseif spectrum == 0

    % 然后，通过 Value 中的情况，人为确定低能态的个数，并获取对应的本征向量
    % 借助 Value 定一个阈值 minus，通过循环将向量取出来
    minus = 1e-63;
    Phi = [];
    for i = 1:L
        if abs(D(i,i)) <= 0.0203
            Phi = [Phi; V(:,i)'];
        end
    end

    % 获得原始坐标对应的态密度 A
    A = zeros(1,L);
    for i = 1:L
        for j = 1:length(Phi(:,1))
            A(i) = Phi(j,i)'*Phi(j,i) + A(i);
        end
    end

    % 将每个原子位置的自由度进行合并，得到该点总的态密度
    S = zeros(1,L/freedom);
    for i = 1:L/freedom
        S(i) = sum(A((i-1)*4+1:(i-1)*4+freedom))/length(Phi(:,1));
    end

    % 接下来就是程序的第二个难点了，就是画图
    % 第一个步骤，是从下到上遍历所有的点，并对每个点获取坐标 (x,y)
    % 我们规定了最左下角的原子的坐标为 (0,0)，并将所有点的坐标存储在 x、y 中
    x = zeros(1,L/freedom);
    y = zeros(1,L/freedom);

    % type 的输入值有 1、2、3、4 四种，对应 AB-C \ C-AB \ AB-AB \ C-C 四种类型
    if type == 1 || type == 3
        if type == 1
            for i = 1:Ny
                x((i-1)*(3*Nx-1)+1:(i-1)*(3*Nx-1)+(2*Nx-1)) = 0:1:(2*Nx-2);
                y((i-1)*(3*Nx-1)+1:(i-1)*(3*Nx-1)+(2*Nx-1)) = 2*(i-1);
                x((i-1)*(3*Nx-1)+2*Nx:i*(3*Nx-1)) = 0:2:(2*Nx-2);
                y((i-1)*(3*Nx-1)+2*Nx:i*(3*Nx-1)) = 2*(i-1) + 1;
            end
            if N_add ~= 0
                x(Ny*(3*Nx-1)+1:L/freedom) = 2*(N_begin-1):1:2*(N_begin+N_add-2);
                y(Ny*(3*Nx-1)+1:L/freedom) = 2*Ny;
            end

        elseif type == 3
            for i = 1:Ny
                x((i-1)*(3*Nx-1)+1:(i-1)*(3*Nx-1)+(2*Nx-1)) = 0:1:(2*Nx-2);
                y((i-1)*(3*Nx-1)+1:(i-1)*(3*Nx-1)+(2*Nx-1)) = 2*(i-1);
                x((i-1)*(3*Nx-1)+2*Nx:i*(3*Nx-1)) = 0:2:(2*Nx-2);
                y((i-1)*(3*Nx-1)+2*Nx:i*(3*Nx-1)) = 2*(i-1) + 1;
            end
            x = x(1:(3*Nx-1)*Ny-Nx);
            y = y(1:(3*Nx-1)*Ny-Nx);
            if N_add ~= 0
                x(Ny*(3*Nx-1)+1-Nx:L/freedom) = 2*(N_begin-1):2:2*(N_begin+N_add-2);
                y(Ny*(3*Nx-1)+1-Nx:L/freedom) = 2*Ny-1;
            end
        end

    elseif type == 2 || type == 4
        if type == 2
            for i = 1:Ny
                x((i-1)*(3*Nx-1)+1:(i-1)*(3*Nx-1)+Nx) = 0:2:(2*Nx-2);
                y((i-1)*(3*Nx-1)+1:(i-1)*(3*Nx-1)+Nx) = 2*(i-1);
                x((i-1)*(3*Nx-1)+Nx+1:i*(3*Nx-1)) = 0:1:(2*Nx-2);
                y((i-1)*(3*Nx-1)+Nx+1:i*(3*Nx-1)) = 2*(i-1) + 1;
            end
            if N_add ~= 0
                x(Ny*(3*Nx-1)+1:L/freedom) = 2*(N_begin-1):2:2*(N_begin+N_add-2);
                y(Ny*(3*Nx-1)+1:L/freedom) = 2*Ny;
            end
            
        elseif type == 4
            for i = 1:Ny
                x((i-1)*(3*Nx-1)+1:(i-1)*(3*Nx-1)+Nx) = 0:2:(2*Nx-2);
                y((i-1)*(3*Nx-1)+1:(i-1)*(3*Nx-1)+Nx) = 2*(i-1);
                x((i-1)*(3*Nx-1)+Nx+1:i*(3*Nx-1)) = 0:1:(2*Nx-2);
                y((i-1)*(3*Nx-1)+Nx+1:i*(3*Nx-1)) = 2*(i-1) + 1;
            end
            x = x(1:(3*Nx-1)*Ny-(2*Nx-1));
            y = y(1:(3*Nx-1)*Ny-(2*Nx-1));
            if N_add ~= 0
                x((3*Nx-1)*Ny-(2*Nx-1)+1:L/freedom) = 2*(N_begin-1):1:2*(N_begin+N_add-2);
                y((3*Nx-1)*Ny-(2*Nx-1)+1:L/freedom) = 2*Ny-1;
            end
        end
    end

    % 根据不同类型的边界条件，画出晶格散点图
    csize = 60;
    linewidth = 1;
    circlewidth = 1;
    scatter(x,y,csize,'k','LineWidth',circlewidth); hold on;
    axis equal;
    axis off;
    set(gcf,'unit','normalized','position',[0.1,0.1,0.8,0.8]);

    % 根据不同的边界条件，将格点间连成线
    if type == 1
        % 先画横线
        for i = 1:Ny
            line([0,2*Nx-2],[2*(i-1),2*(i-1)],'Color','black','LineWidth',linewidth); hold on;
        end
        % 再画纵线
        for i = 1:Nx
            line([(i-1)*2,(i-1)*2],[0,2*Ny-1],'Color','black','LineWidth',linewidth); hold on;
        end
        if N_add ~= 0
            line([2*(N_begin-1),2*(N_begin+N_add-2)],[2*Ny,2*Ny],'Color','black','LineWidth',linewidth); hold on;
            for j = 1:N_add
                line([2*(N_begin-1)+2*(j-1),2*(N_begin-1)+2*(j-1)],[2*Ny-1,2*Ny],'Color','black','LineWidth',linewidth); hold on;
            end
        end
    elseif type == 2
        % 先画横线
        for i = 1:Ny
            line([0,2*Nx-2],[2*i-1,2*i-1],'Color','black','LineWidth',linewidth); hold on;
        end
        % 再画纵线
        for i = 1:Nx
            line([(i-1)*2,(i-1)*2],[0,2*Ny-1],'Color','black','LineWidth',linewidth); hold on;
        end
        if N_add ~= 0
            for j = 1:N_add
                line([2*(N_begin-1)+2*(j-1),2*(N_begin-1)+2*(j-1)],[2*Ny-1,2*Ny],'Color','black','LineWidth',linewidth); hold on;
            end
        end
    elseif type == 3
        % 先画横线
        for i = 1:Ny
            line([0,2*Nx-2],[2*(i-1),2*(i-1)],'Color','black','LineWidth',linewidth); hold on;
        end
        % 再画纵线
        for i = 1:Nx
            line([(i-1)*2,(i-1)*2],[0,2*Ny-2],'Color','black','LineWidth',linewidth); hold on;
        end
        if N_add ~= 0
            for j = 1:N_add
                line([2*(N_begin-1)+2*(j-1),2*(N_begin-1)+2*(j-1)],[2*Ny-2,2*Ny-1],'Color','black','LineWidth',linewidth); hold on;
            end
        end
    elseif type == 4
        % 先画横线
        for i = 1:Ny-1
            line([0,2*Nx-2],[2*i-1,2*i-1],'Color','black','LineWidth',linewidth); hold on;
        end
        % 再画纵线
        for i = 1:Nx
            line([(i-1)*2,(i-1)*2],[0,2*Ny-2],'Color','black','LineWidth',linewidth); hold on;
        end
        if N_add ~= 0
            line([2*(N_begin-1),2*(N_begin+N_add-2)],[2*Ny-1,2*Ny-1],'Color','black','LineWidth',linewidth); hold on;
            for j = 1:N_add
                line([2*(N_begin-1)+2*(j-1),2*(N_begin-1)+2*(j-1)],[2*Ny-2,2*Ny-1],'Color','black','LineWidth',linewidth); hold on;
            end
        end
    end

    % 接下来，给对应原子涂色
    % 这里分两种画图，一种是示意图，一种是态密度图
    % lattice = 1 代表示意图；lattice = 0 代表态密度图
    scatter(x,y,csize,'white','filled');

    if lattice == 0
        Max = max(S);
        C(:,2) = 1/Max*(Max-S(:));
        C(:,3) = 1/Max*(Max-S(:));
        C(:,1) = 1; % 红色分布
        scatter(x,y,csize,C,'filled');
    elseif lattice ~= 0
        if type == 1 || type == 3
            for i = 1:length(x)
                % A 类原子着色
                if mod(x(i),2) == 0 && mod(y(i),2) == 0
                    scatter(x(i),y(i),csize,'r','filled');
                    % B 类原子着色
                elseif mod(x(i),2) == 1 && mod(y(i),2) == 0
                    scatter(x(i),y(i),csize,[0.5,0,0.7],'filled');
                    % C 类原子着色
                elseif mod(x(i),2) == 0 && mod(y(i),2) == 1
                    scatter(x(i),y(i),csize,'b','filled');
                end
            end
        elseif type == 2 || type == 4
            for i = 1:length(x)
                % A 类原子着色
                if mod(x(i),2) == 0 && mod(y(i),2) == 1
                    scatter(x(i),y(i),csize,'r','filled');
                    % B 类原子着色
                elseif mod(x(i),2) == 1 && mod(y(i),2) == 1
                    scatter(x(i),y(i),csize,[0.5,0,0.7],'filled');
                    % C 类原子着色
                elseif mod(x(i),2) == 0 && mod(y(i),2) == 0
                    scatter(x(i),y(i),csize,'b','filled');
                end
            end
        end
    end
    % set(gcf,'unit','normalized','position',[0.1,0.1,0.8,0.8]);
    saveas(gcf,str2);

end


