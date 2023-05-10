function [Energy] = Energyspectrum(k)

% 本程序用于计算已知哈密顿量的能谱图

% 首先导出哈密顿量并对角化
H = Hamiltonian(k) ;
[V,D] = eig(H);                         % V 是本征向量，D 是本征值
[Y,I] = sort(diag(D),'descend');        % 降序，Y 是已经排序好的数组，而 I 是 Y 每个数在原来数组（矩阵）中的位置

% 直接将 Y 作为一个列向量输出
Energy = Y ;

end