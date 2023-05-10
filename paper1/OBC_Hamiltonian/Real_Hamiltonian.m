function[H] = Real_Hamiltonian(Nx,Ny,type)

% type 代表边界类型，type = 1 是第一类，type = 2 是第二类
% Nx、Ny 代表沿 x 和 y 方向具有的初基原胞 A-B

t = 1;
lambda_SO = 0.1;
lambda_v  = 0.1;
gamma = 0.4;
Mx = 0.0;
My = 0.0;

sigma0 = [ 1,  0;   0,  1];
sigmax = [ 0,  1;   1,  0];
sigmay = [ 0, 1i; -1i,  0];
sigmaz = [ 1,  0;   0, -1];

judge = 0;
if mod(Ny,2) == 0
    Ny = Ny / 2;
elseif mod(Ny,2) == 1
    Ny = (Ny+1)/2;
    judge = 1;
end

if type == 2
    Ny = Ny + 1;
end

% 扩大晶胞，为了保证索引为正
H = zeros((Nx+2)*(Ny+2)*2*2,(Nx+2)*(Ny+2)*2*2);

% H 中的标记分别是：
% A11, A12, ... A1Nx, B11, B12, ... B1Nx, ... ANy1, ANy2, ... ANyNx, BNy1, BNy2, ... BNyNx
% 每个坐标还会存储上自旋和下自旋两种态，因此维度会翻倍
% 不考虑自旋效应下
% A(x,y) 在 H 中的坐标为 ( x+Nx*(2*y-2) )
% B(x,y) 在 H 中的坐标为 ( x+Nx*(2*y-1) )
% 考虑自旋效应，还得将坐标扩展一倍
% 注意：索引可能不存在（小于 0 ），这个要额外处理一下
% 采用的思路是扩展晶胞一圈，最后再去除

% 基础作用项
H0 = zeros((Nx+2)*(Ny+2)*2,(Nx+2)*(Ny+2)*2);

for x = 2:Nx+1
    for y = 2:Ny+1
        if mod(y,2) == 1
            % 以 A 原子为例，研究 A-B 相互作用，B 只需取共轭
            H0(x+(Nx+2)*(2*y-2), x-1+(Nx+2)*(2*y-1)) = t;
            H0(x+(Nx+2)*(2*y-2), x  +(Nx+2)*(2*y-1)) = t;
            H0(x+(Nx+2)*(2*y-2), x  +(Nx+2)*(2*y-3)) = t;
        else
            H0(x+(Nx+2)*(2*y-2), x  +(Nx+2)*(2*y-1)) = t;
            H0(x+(Nx+2)*(2*y-2), x+1+(Nx+2)*(2*y-1)) = t;
            H0(x+(Nx+2)*(2*y-2), x  +(Nx+2)*(2*y-3)) = t;
        end
    end
end
H0 = H0 + H0';
% 在考虑自旋自由度时，矩阵的维度会翻倍
H0 = kron(H0, sigma0);


% SOC 相互作用
H_SO = zeros((Nx+2)*(Ny+2)*2,(Nx+2)*(Ny+2)*2);

for x = 2:Nx+1
    for y = 2:Ny+1
        if mod(y,2) == 1
            % 对 A 的次近邻相互作用的分析
            H_SO(x+(Nx+2)*(2*y-2), x+1+(Nx+2)*(2*y-2)) = -lambda_SO; % 右
            H_SO(x+(Nx+2)*(2*y-2), x-1+(Nx+2)*(2*y-2)) =  lambda_SO; % 左
            H_SO(x+(Nx+2)*(2*y-2), x  +(Nx+2)*(2*y  )) =  lambda_SO; % 右上
            H_SO(x+(Nx+2)*(2*y-2), x-1+(Nx+2)*(2*y  )) = -lambda_SO; % 左上
            H_SO(x+(Nx+2)*(2*y-2), x  +(Nx+2)*(2*y-4)) =  lambda_SO; % 右下
            H_SO(x+(Nx+2)*(2*y-2), x-1+(Nx+2)*(2*y-4)) = -lambda_SO; % 左下
            % 对 B 的次近邻相互作用的分析
            H_SO(x+(Nx+2)*(2*y-1), x+1+(Nx+2)*(2*y-1)) =  lambda_SO; % 右
            H_SO(x+(Nx+2)*(2*y-1), x-1+(Nx+2)*(2*y-1)) = -lambda_SO; % 左
            H_SO(x+(Nx+2)*(2*y-1), x+1+(Nx+2)*(2*y+1)) = -lambda_SO; % 右上
            H_SO(x+(Nx+2)*(2*y-1), x  +(Nx+2)*(2*y+1)) =  lambda_SO; % 左上
            H_SO(x+(Nx+2)*(2*y-1), x+1+(Nx+2)*(2*y-3)) = -lambda_SO; % 右下
            H_SO(x+(Nx+2)*(2*y-1), x  +(Nx+2)*(2*y-3)) =  lambda_SO; % 左下
        elseif mod(y,2) == 0
            % 对 A 的次近邻相互作用的分析
            H_SO(x+(Nx+2)*(2*y-2), x+1+(Nx+2)*(2*y-2)) = -lambda_SO; % 右
            H_SO(x+(Nx+2)*(2*y-2), x-1+(Nx+2)*(2*y-2)) =  lambda_SO; % 左
            H_SO(x+(Nx+2)*(2*y-2), x+1+(Nx+2)*(2*y  )) =  lambda_SO; % 右上
            H_SO(x+(Nx+2)*(2*y-2), x  +(Nx+2)*(2*y  )) = -lambda_SO; % 左上
            H_SO(x+(Nx+2)*(2*y-2), x+1+(Nx+2)*(2*y-4)) =  lambda_SO; % 右下
            H_SO(x+(Nx+2)*(2*y-2), x  +(Nx+2)*(2*y-4)) = -lambda_SO; % 左下
            % 对 B 的次近邻相互作用的分析
            H_SO(x+(Nx+2)*(2*y-1), x+1+(Nx+2)*(2*y-1)) =  lambda_SO; % 右
            H_SO(x+(Nx+2)*(2*y-1), x-1+(Nx+2)*(2*y-1)) = -lambda_SO; % 左
            H_SO(x+(Nx+2)*(2*y-1), x  +(Nx+2)*(2*y+1)) = -lambda_SO; % 右上
            H_SO(x+(Nx+2)*(2*y-1), x-1+(Nx+2)*(2*y+1)) =  lambda_SO; % 左上
            H_SO(x+(Nx+2)*(2*y-1), x  +(Nx+2)*(2*y-3)) = -lambda_SO; % 右下
            H_SO(x+(Nx+2)*(2*y-1), x-1+(Nx+2)*(2*y-3)) =  lambda_SO; % 左下
        end
    end
end

H_SO = 1i*H_SO;

% 由于上下自旋对应的态不一样，差一个负号
H_SO = kron(H_SO,sigmaz);


% 化学势项

H_v = zeros((Nx+2)*(Ny+2)*2,(Nx+2)*(Ny+2)*2);

for x = 1:Nx+2
    for y = 1:Ny+2
        % 对 A 的势能项
        H_v(x+(Nx+2)*(2*y-2),x+(Nx+2)*(2*y-2)) =  lambda_v;
        % 对 B 的势能项
        H_v(x+(Nx+2)*(2*y-1),x+(Nx+2)*(2*y-1)) = -lambda_v;
    end
end

% 在考虑自旋自由度时，矩阵的维度会翻倍
H_v = kron(H_v,sigma0);


% 磁场相互作用
H_M = zeros((Nx+2)*(Ny+2)*2,(Nx+2)*(Ny+2)*2);

for x = 1:Nx+2
    for y = 1:Ny+2
        % 对 A 的势能项
        H_M(x+(Nx+2)*(2*y-2),x+(Nx+2)*(2*y-2)) = 1;
        % 对 B 的势能项
        H_M(x+(Nx+2)*(2*y-1),x+(Nx+2)*(2*y-1)) = gamma;
    end
end
M = Mx*sigmax + My*sigmay;
H_M = kron(H_M, M);


% 总的哈密顿量为
H = H0 + H_SO + H_v + H_M;

% 由于采用了扩大晶格的办法，因此需要将边缘一圈去除掉，取出中间的那部分晶格
% 先去除上面和下面的两层，取出中间夹层部分
H = H(((Nx+2)*4+1):((Nx+2)*(Ny+1)*4),((Nx+2)*4+1):((Nx+2)*(Ny+1)*4));
% 然后再将左右两边的边缘去除，此时有 Ny 层，(Nx+2) 列
% 为了避免索引出错，应该从坐标量大的开始去除
% A(x,y) 在 H 中的坐标为 ( x+(Nx+2)*(2*y-2) )
% B(x,y) 在 H 中的坐标为 ( x+(Nx+2)*(2*y-1) )

for s = 1:Ny
    y = Ny - s + 1;
    H(((Nx+2)*2-1+(Nx+2)*(2*y-1)*2):((Nx+2)*2+(Nx+2)*(2*y-1)*2),:) = [];
    H(:,((Nx+2)*2-1+(Nx+2)*(2*y-1)*2):((Nx+2)*2+(Nx+2)*(2*y-1)*2)) = [];
    H((1+(Nx+2)*(2*y-1)*2):(2+(Nx+2)*(2*y-1)*2),:) = [];
    H(:,(1+(Nx+2)*(2*y-1)*2):(2+(Nx+2)*(2*y-1)*2)) = [];
    H(((Nx+2)*2-1+(Nx+2)*(2*y-2)*2):((Nx+2)*2+(Nx+2)*(2*y-2)*2),:) = [];
    H(:,((Nx+2)*2-1+(Nx+2)*(2*y-2)*2):((Nx+2)*2+(Nx+2)*(2*y-2)*2)) = [];
    H((1+(Nx+2)*(2*y-2)*2):(2+(Nx+2)*(2*y-2)*2),:) = [];
    H(:,(1+(Nx+2)*(2*y-2)*2):(2+(Nx+2)*(2*y-2)*2)) = [];
end

% 接下来，对边界类型进行分类
% type = 1 对应的是边界均为 A-A/B 型，哈密顿量不需要改变
% type = 2 对应的是边界均为 B-A/B 型，哈密顿量需要去掉最下面的 A 层

if type == 2
    H((1+2*Nx*(2*Ny-1)):(4*Nx*Ny),:) = [];
    H(:,(1+2*Nx*(2*Ny-1)):(4*Nx*Ny)) = [];
    H(1:Nx*2,:) = [];
    H(:,1:Nx*2) = [];
end

if judge == 1
    L = length(H(1,:));
    H(L-Nx*2+1:L,:) = [];
    H(:,L-Nx*2+1:L) = [];
end

% % 此时 type = 1 得到的 Matrix 是 Nx*Ny*2*2 维度
% % 此时 type = 2 得到的 Matrix 是 Nx*(Ny-1)*2*2 维度
% % 接下来，为了让系统左右对称，我们需要去除右边的一部分
% % 具体去除的元素为：
% % I  类边界是去除右侧的 B1Nx, A2Nx, B3Nx, A4Nx, ... (注意，第一行是A)
% % II 类边界是去除左侧的 B21, A31, B41, A51, ...(注意，第一行是B)
%
% if type == 1
%     for y = 1:2*Ny
%         s = Ny + 1 - y
%         if mod(s,4) == 0
%             H(2*Nx*(s-1+1),:) = [];
%             H(2*Nx*(s-2+1),:) = [];
%         end
%     end
% elseif type == 2
%     for y = 1:2*Ny-2
%         s = Ny + 1 - y
%         if mod(s,4) == 0
%             H(2*Nx*(s)+1,:) = [];
%             H(2*Nx*(s-1)+1,:) = [];
%         end
%     end
% end

end