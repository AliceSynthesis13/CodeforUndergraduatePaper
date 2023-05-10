function [Energy] = Energyspectrum(k,n,type)

%本程序用于计算已知哈密顿量的能谱图

%首先导出哈密顿量并对角化
H = OBC_Hamiltonian(k,n,type) ;
[V,D] = eig(H);                         %V是本征向量，D是本征值
[Y,I] = sort(diag(D),'descend');        %降序，Y是已经排序好的数组，而I是Y每个数在原来数组（矩阵）中的位置

%直接讲Y作为一个列向量输出
Energy = Y ;

end