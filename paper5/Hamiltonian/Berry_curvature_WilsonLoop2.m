function [Berry_curvature_WilsonLoop2] = Berry_curvature_WilsonLoop2(R,dkx,dky,rows)
%% 解释面板
% 本程序用于计算所有哈密顿量为矩阵的 Berry Curvature , 并且适用于能谱简并的哈密顿量
% 接下来阐述几个输入参数的意义：
%   R:输入一个任意维度的参数向量，使用者必须记录清楚其意义

% 本程序选择围绕定点 R 周围的 8 个点的 Loop 进行计算，因为我们认为dkx和dky是足够小以至于可以忽略的【微分量】，
% 那么环路积分理所应当的转化成这八个点的乘积


%% 第一步 计算围绕 R 点的 8 个点的特征向量 ， 从左下角开始逆时针编号 1 ～ 8 

[V1,~] = eig(Hamiltonian(R+[   0,   0]));
[V2,~] = eig(Hamiltonian(R+[+dkx,   0]));
[V3,~] = eig(Hamiltonian(R+[+dkx,+dky]));
[V4,~] = eig(Hamiltonian(R+[   0,+dky]));



%% 第二步 计算 Loop 的内积乘积
% 此处特注：我们计算的是所有能带的结果，但只有对角元对应着 Berry phase

phi = - imag( log( (V1'*V2.*eye(rows,rows)) .* (V2'*V3.*eye(rows,rows)) .* (V3'*V4.*eye(rows,rows)) .* (V4'*V1.*eye(rows,rows)) ) );

Berry_curvature_WilsonLoop2 = diag(phi/(dkx*dky));



end
