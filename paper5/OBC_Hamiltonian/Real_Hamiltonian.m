function[H] = Real_Hamiltonian(Nx,Ny,type,t,lambda_SO,mu,Delta)

% t = 1;
% lambda_SO = 0.1;
% mu = 0.0;
% Delta = 0.0;

% 这玩意是实空间的哈密顿量，输入晶格相关参数 Nx、Ny、type
% type 的输入值有 1、2、3、4 四种，对应 AB-C \ C-AB \ AB-AB \ C-C 四种类型
% 为了避免边界处索引出错，这里还是采用扩展晶胞再删除的办法
Nx = Nx + 2;
Ny = Ny + 2;

sigma0 = [ 1,  0;   0,  1];
sigmax = [ 0,  1;   1,  0];
sigmay = [ 0, 1i; -1i,  0];
sigmaz = [ 1,  0;   0, -1];

freedom = 4; % 代表每个格点具有的自由度

H = zeros(Ny*Nx*3*freedom,Ny*Nx*3*freedom);
% 在 j = Ny 处的子晶格自由度可能有 AB 和 C 两种类型，这里采用删除法来进行变换
% 对于 (1)AB-C \ (3)AB-AB 只需要在 (1) 的基础上删去最上面的 C  原子
% 对于 (2)C-AB \ (4)C-C   只需要在 (2) 的基础上删去最上面的 AB 原子

% posi 函数的输入形式是：posi(i,j,sublattice,type,Nx)

% 接下来计算第一项相互作用（近邻相互作用）
H0 = zeros(Ny*Nx*3*freedom,Ny*Nx*3*freedom);
A = kron(sigmaz,sigma0);
% A = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0];
if type == 1 || type == 3
    for i = 2:(Nx-1)
        for j = 2:(Ny-1)
            % A-B
            H0(posi(i,j,1,type,Nx)+1:posi(i,j,1,type,Nx)+4, posi(i-1,j,2,type,Nx)+1:posi(i-1,j,2,type,Nx)+4) = -t*A;
            H0(posi(i,j,1,type,Nx)+1:posi(i,j,1,type,Nx)+4, posi(i  ,j,2,type,Nx)+1:posi(i  ,j,2,type,Nx)+4) = -t*A;
            % B-A
            H0(posi(i,j,2,type,Nx)+1:posi(i,j,2,type,Nx)+4, posi(i  ,j,1,type,Nx)+1:posi(i  ,j,1,type,Nx)+4) = -t*A;
            H0(posi(i,j,2,type,Nx)+1:posi(i,j,2,type,Nx)+4, posi(i+1,j,1,type,Nx)+1:posi(i+1,j,1,type,Nx)+4) = -t*A;
            % A-C
            H0(posi(i,j,1,type,Nx)+1:posi(i,j,1,type,Nx)+4, posi(i,j-1,3,type,Nx)+1:posi(i,j-1,3,type,Nx)+4) = -t*A;
            H0(posi(i,j,1,type,Nx)+1:posi(i,j,1,type,Nx)+4, posi(i,j  ,3,type,Nx)+1:posi(i,j  ,3,type,Nx)+4) = -t*A;
            % C-A
            H0(posi(i,j,3,type,Nx)+1:posi(i,j,3,type,Nx)+4, posi(i,j  ,1,type,Nx)+1:posi(i,j  ,1,type,Nx)+4) = -t*A;
            H0(posi(i,j,3,type,Nx)+1:posi(i,j,3,type,Nx)+4, posi(i,j+1,1,type,Nx)+1:posi(i,j+1,1,type,Nx)+4) = -t*A;
        end
    end
elseif type == 2 || type == 4
    for i = 2:(Nx-1)
        for j = 2:(Ny-1)
            % A-B
            H0(posi(i,j,1,type,Nx)+1:posi(i,j,1,type,Nx)+4, posi(i-1,j,2,type,Nx)+1:posi(i-1,j,2,type,Nx)+4) = -t*A;
            H0(posi(i,j,1,type,Nx)+1:posi(i,j,1,type,Nx)+4, posi(i  ,j,2,type,Nx)+1:posi(i  ,j,2,type,Nx)+4) = -t*A;
            % B-A
            H0(posi(i,j,2,type,Nx)+1:posi(i,j,2,type,Nx)+4, posi(i  ,j,1,type,Nx)+1:posi(i  ,j,1,type,Nx)+4) = -t*A;
            H0(posi(i,j,2,type,Nx)+1:posi(i,j,2,type,Nx)+4, posi(i+1,j,1,type,Nx)+1:posi(i+1,j,1,type,Nx)+4) = -t*A;
            % A-C
            H0(posi(i,j,1,type,Nx)+1:posi(i,j,1,type,Nx)+4, posi(i,j  ,3,type,Nx)+1:posi(i,j  ,3,type,Nx)+4) = -t*A;
            H0(posi(i,j,1,type,Nx)+1:posi(i,j,1,type,Nx)+4, posi(i,j+1,3,type,Nx)+1:posi(i,j+1,3,type,Nx)+4) = -t*A;
            % C-A
            H0(posi(i,j,3,type,Nx)+1:posi(i,j,3,type,Nx)+4, posi(i,j-1,1,type,Nx)+1:posi(i,j-1,1,type,Nx)+4) = -t*A;
            H0(posi(i,j,3,type,Nx)+1:posi(i,j,3,type,Nx)+4, posi(i,j  ,1,type,Nx)+1:posi(i,j  ,1,type,Nx)+4) = -t*A;
        end
    end
end

% 接下来算第二项（次近邻相互作用，SOC 相互作用）
H_SO = zeros(Ny*Nx*3*freedom,Ny*Nx*3*freedom);
B = kron(sigma0,sigmaz);
if type == 1 || type == 3
    for i = 2:(Nx-1)
        for j = 2:(Ny-1)
            % B to C
            H_SO(posi(i,j,2,type,Nx)+1:posi(i,j,2,type,Nx)+4, posi(i  ,j  ,3,type,Nx)+1:posi(i  ,j  ,3,type,Nx)+4) = -1i*lambda_SO*B;
            H_SO(posi(i,j,2,type,Nx)+1:posi(i,j,2,type,Nx)+4, posi(i+1,j  ,3,type,Nx)+1:posi(i+1,j  ,3,type,Nx)+4) =  1i*lambda_SO*B;
            H_SO(posi(i,j,2,type,Nx)+1:posi(i,j,2,type,Nx)+4, posi(i  ,j-1,3,type,Nx)+1:posi(i  ,j-1,3,type,Nx)+4) =  1i*lambda_SO*B;
            H_SO(posi(i,j,2,type,Nx)+1:posi(i,j,2,type,Nx)+4, posi(i+1,j-1,3,type,Nx)+1:posi(i+1,j-1,3,type,Nx)+4) = -1i*lambda_SO*B;
            % C to B
            H_SO(posi(i,j,3,type,Nx)+1:posi(i,j,3,type,Nx)+4, posi(i-1,j+1,2,type,Nx)+1:posi(i-1,j+1,2,type,Nx)+4) =  1i*lambda_SO*B;
            H_SO(posi(i,j,3,type,Nx)+1:posi(i,j,3,type,Nx)+4, posi(i  ,j+1,2,type,Nx)+1:posi(i  ,j+1,2,type,Nx)+4) = -1i*lambda_SO*B;
            H_SO(posi(i,j,3,type,Nx)+1:posi(i,j,3,type,Nx)+4, posi(i-1,j  ,2,type,Nx)+1:posi(i-1,j  ,2,type,Nx)+4) = -1i*lambda_SO*B;
            H_SO(posi(i,j,3,type,Nx)+1:posi(i,j,3,type,Nx)+4, posi(i  ,j  ,2,type,Nx)+1:posi(i  ,j  ,2,type,Nx)+4) =  1i*lambda_SO*B;
        end
    end
elseif type == 2 || type == 4
    for i = 2:(Nx-1)
        for j = 2:(Ny-1)
            % B to C
            H_SO(posi(i,j,2,type,Nx)+1:posi(i,j,2,type,Nx)+4, posi(i  ,j+1,3,type,Nx)+1:posi(i  ,j+1,3,type,Nx)+4) = -1i*lambda_SO*B;
            H_SO(posi(i,j,2,type,Nx)+1:posi(i,j,2,type,Nx)+4, posi(i+1,j+1,3,type,Nx)+1:posi(i+1,j+1,3,type,Nx)+4) =  1i*lambda_SO*B;
            H_SO(posi(i,j,2,type,Nx)+1:posi(i,j,2,type,Nx)+4, posi(i  ,j  ,3,type,Nx)+1:posi(i  ,j  ,3,type,Nx)+4) =  1i*lambda_SO*B;
            H_SO(posi(i,j,2,type,Nx)+1:posi(i,j,2,type,Nx)+4, posi(i+1,j  ,3,type,Nx)+1:posi(i+1,j  ,3,type,Nx)+4) = -1i*lambda_SO*B;
            % C to B
            H_SO(posi(i,j,3,type,Nx)+1:posi(i,j,3,type,Nx)+4, posi(i-1,j  ,2,type,Nx)+1:posi(i-1,j  ,2,type,Nx)+4) =  1i*lambda_SO*B;
            H_SO(posi(i,j,3,type,Nx)+1:posi(i,j,3,type,Nx)+4, posi(i  ,j  ,2,type,Nx)+1:posi(i  ,j  ,2,type,Nx)+4) = -1i*lambda_SO*B;
            H_SO(posi(i,j,3,type,Nx)+1:posi(i,j,3,type,Nx)+4, posi(i-1,j-1,2,type,Nx)+1:posi(i-1,j-1,2,type,Nx)+4) = -1i*lambda_SO*B;
            H_SO(posi(i,j,3,type,Nx)+1:posi(i,j,3,type,Nx)+4, posi(i  ,j-1,2,type,Nx)+1:posi(i  ,j-1,2,type,Nx)+4) =  1i*lambda_SO*B;
        end
    end
end


% 接下来算第三项，化学势
H_mu = zeros(Ny*Nx*3*freedom,Ny*Nx*3*freedom);
C = kron(sigmaz,sigma0);
for i = 2:(Nx-1)
    for j = 2:(Ny-1)
        H_mu(posi(i,j,1,type,Nx)+1:posi(i,j,1,type,Nx)+4,posi(i,j,1,type,Nx)+1:posi(i,j,1,type,Nx)+4) = -mu*C;
        H_mu(posi(i,j,2,type,Nx)+1:posi(i,j,2,type,Nx)+4,posi(i,j,2,type,Nx)+1:posi(i,j,2,type,Nx)+4) = -mu*C;
        H_mu(posi(i,j,3,type,Nx)+1:posi(i,j,3,type,Nx)+4,posi(i,j,3,type,Nx)+1:posi(i,j,3,type,Nx)+4) = -mu*C;
    end
end


% 接下来计算第四项，d 波的超导配对
H_SC = zeros(Ny*Nx*3*freedom,Ny*Nx*3*freedom);
D = kron(sigmay,sigmay);
if type == 1 || type == 3
    for i = 2:(Nx-1)
        for j = 2:(Ny-1)
            % A-B
            H_SC(posi(i,j,1,type,Nx)+1:posi(i,j,1,type,Nx)+4, posi(i-1,j,2,type,Nx)+1:posi(i-1,j,2,type,Nx)+4) =  Delta*D;
            H_SC(posi(i,j,1,type,Nx)+1:posi(i,j,1,type,Nx)+4, posi(i  ,j,2,type,Nx)+1:posi(i  ,j,2,type,Nx)+4) =  Delta*D;
            % B-A
            H_SC(posi(i,j,2,type,Nx)+1:posi(i,j,2,type,Nx)+4, posi(i  ,j,1,type,Nx)+1:posi(i  ,j,1,type,Nx)+4) =  Delta*D;
            H_SC(posi(i,j,2,type,Nx)+1:posi(i,j,2,type,Nx)+4, posi(i+1,j,1,type,Nx)+1:posi(i+1,j,1,type,Nx)+4) =  Delta*D;
            % A-C
            H_SC(posi(i,j,1,type,Nx)+1:posi(i,j,1,type,Nx)+4, posi(i,j-1,3,type,Nx)+1:posi(i,j-1,3,type,Nx)+4) = -Delta*D;
            H_SC(posi(i,j,1,type,Nx)+1:posi(i,j,1,type,Nx)+4, posi(i,j  ,3,type,Nx)+1:posi(i,j  ,3,type,Nx)+4) = -Delta*D;
            % C-A
            H_SC(posi(i,j,3,type,Nx)+1:posi(i,j,3,type,Nx)+4, posi(i,j  ,1,type,Nx)+1:posi(i,j  ,1,type,Nx)+4) = -Delta*D;
            H_SC(posi(i,j,3,type,Nx)+1:posi(i,j,3,type,Nx)+4, posi(i,j+1,1,type,Nx)+1:posi(i,j+1,1,type,Nx)+4) = -Delta*D;
        end
    end
elseif type == 2 || type == 4
    for i = 2:(Nx-1)
        for j = 2:(Ny-1)
            % A-B
            H_SC(posi(i,j,1,type,Nx)+1:posi(i,j,1,type,Nx)+4, posi(i-1,j,2,type,Nx)+1:posi(i-1,j,2,type,Nx)+4) =  Delta*D;
            H_SC(posi(i,j,1,type,Nx)+1:posi(i,j,1,type,Nx)+4, posi(i  ,j,2,type,Nx)+1:posi(i  ,j,2,type,Nx)+4) =  Delta*D;
            % B-A
            H_SC(posi(i,j,2,type,Nx)+1:posi(i,j,2,type,Nx)+4, posi(i  ,j,1,type,Nx)+1:posi(i  ,j,1,type,Nx)+4) =  Delta*D;
            H_SC(posi(i,j,2,type,Nx)+1:posi(i,j,2,type,Nx)+4, posi(i+1,j,1,type,Nx)+1:posi(i+1,j,1,type,Nx)+4) =  Delta*D;
            % A-C
            H_SC(posi(i,j,1,type,Nx)+1:posi(i,j,1,type,Nx)+4, posi(i,j  ,3,type,Nx)+1:posi(i,j  ,3,type,Nx)+4) = -Delta*D;
            H_SC(posi(i,j,1,type,Nx)+1:posi(i,j,1,type,Nx)+4, posi(i,j+1,3,type,Nx)+1:posi(i,j+1,3,type,Nx)+4) = -Delta*D;
            % C-A
            H_SC(posi(i,j,3,type,Nx)+1:posi(i,j,3,type,Nx)+4, posi(i,j-1,1,type,Nx)+1:posi(i,j-1,1,type,Nx)+4) = -Delta*D;
            H_SC(posi(i,j,3,type,Nx)+1:posi(i,j,3,type,Nx)+4, posi(i,j  ,1,type,Nx)+1:posi(i,j  ,1,type,Nx)+4) = -Delta*D;
        end
    end
end

% 总的哈密顿量为：
H = (1/2)*(H0 + H_SO + H_mu + H_SC);
% H = H0 + H_SC + H_mu;

% 接下来将外面一圈去除，这样保证相互作用不会少项数的情况下得到正确的哈密顿量
% 为了保证顺序不会出错，我们应该从索引值高的往低的进行
% 先去除上面一层和下面一层
if type == 1 || type == 3
    H(posi(1,Ny,1,type,Nx)+1:posi(Nx,Ny,3,type,Nx)+4,:) = []; % 上边界
    H(:,posi(1,Ny,1,type,Nx)+1:posi(Nx,Ny,3,type,Nx)+4) = []; % 上边界
    H(posi(1, 1,1,type,Nx)+1:posi(Nx, 1,3,type,Nx)+4,:) = []; % 下边界
    H(:,posi(1, 1,1,type,Nx)+1:posi(Nx, 1,3,type,Nx)+4) = []; % 下边界
elseif type == 2 || type == 4
    H(posi(1,Ny,3,type,Nx)+1:posi(Nx,Ny,2,type,Nx)+4,:) = []; % 上边界
    H(:,posi(1,Ny,3,type,Nx)+1:posi(Nx,Ny,2,type,Nx)+4) = []; % 上边界
    H(posi(1, 1,3,type,Nx)+1:posi(Nx, 1,2,type,Nx)+4,:) = []; % 下边界
    H(:,posi(1, 1,3,type,Nx)+1:posi(Nx, 1,2,type,Nx)+4) = []; % 下边界
end

% 接下来去除左右两边的边界
% 左边需要去掉 ABC，右边需要去掉 B+ABC

Ny = Ny - 2;
if type == 1 || type == 3
    for j = 1:Ny
        y = Ny + 1 - j;
        H(posi(Nx,y,3,type,Nx)+1:posi(Nx,y,3,type,Nx)+4,:) = [];
        H(:,posi(Nx,y,3,type,Nx)+1:posi(Nx,y,3,type,Nx)+4) = [];
        H(posi(1,y,3,type,Nx)+1:posi(1,y,3,type,Nx)+4,:) = [];
        H(:,posi(1,y,3,type,Nx)+1:posi(1,y,3,type,Nx)+4) = [];
        H(posi(Nx-1,y,2,type,Nx)+1:posi(Nx,y,2,type,Nx)+4,:) = [];
        H(:,posi(Nx-1,y,2,type,Nx)+1:posi(Nx,y,2,type,Nx)+4) = [];
        H(posi(1,y,1,type,Nx)+1:posi(1,y,2,type,Nx)+4,:) = [];
        H(:,posi(1,y,1,type,Nx)+1:posi(1,y,2,type,Nx)+4) = [];
    end
elseif type == 2 || type == 4
    for j = 1:Ny
        y = Ny + 1 - j;
        H(posi(Nx-1,y,2,type,Nx)+1:posi(Nx,y,2,type,Nx)+4,:) = [];
        H(:,posi(Nx-1,y,2,type,Nx)+1:posi(Nx,y,2,type,Nx)+4) = [];
        H(posi(1,y,1,type,Nx)+1:posi(1,y,2,type,Nx)+4,:) = [];
        H(:,posi(1,y,1,type,Nx)+1:posi(1,y,2,type,Nx)+4) = [];
        H(posi(Nx,y,3,type,Nx)+1:posi(Nx,y,3,type,Nx)+4,:) = [];
        H(:,posi(Nx,y,3,type,Nx)+1:posi(Nx,y,3,type,Nx)+4) = [];
        H(posi(1,y,3,type,Nx)+1:posi(1,y,3,type,Nx)+4,:) = [];
        H(:,posi(1,y,3,type,Nx)+1:posi(1,y,3,type,Nx)+4) = [];       
    end
end

Nx = Nx - 2;
% 这个时候，四种哈密顿量的维度是相同的，我们接下来需要针对不同类型进行处理
% 注意此时的 Nx 和 Ny 又变得和输入的值相同了，x 方向有 Nx 个 A
% 但是，由于我们将边界处理了，左右两边不同，因此不能用原先的 posi 函数来计算了
% 重新定义了函数 posim
% type 的输入值有 1、2、3、4 四种，对应 AB-C \ C-AB \ AB-AB \ C-C 四种类型
% 对 type = 1 和 type = 2 这两种，是不需要做任何变更的
% 对 type = 3 和 type = 4 这两种，需要删除最上面一行
% 对 type = 3 删除的是最上面的 C  原子层
% 对 type = 4 删除的是最上面的 AB 原子层

if type == 3
    H(posim(1,Ny,3,type,Nx)+1:posim(Nx,Ny,3,type,Nx)+4,:) = [];
    H(:,posim(1,Ny,3,type,Nx)+1:posim(Nx,Ny,3,type,Nx)+4) = [];
elseif type == 4
    H(posim(1,Ny,1,type,Nx)+1:posim(Nx,Ny,1,type,Nx)+4,:) = [];
    H(:,posim(1,Ny,1,type,Nx)+1:posim(Nx,Ny,1,type,Nx)+4) = [];
end


end


