function[s] = posim(i,j,sublattice,type,Nx)

% 这个函数实现的功能是：输入初基原胞的坐标和子晶格自由度，返回矩阵对应元的位置
% i 的取值是从 1 到 Nx，j 的取值是从 1 到 Ny
% Nx 方向边界保持为 AC-AC 型，则 Nx 代表 A 类或 C 类原子的横向个数
% Ny 方向边界是有三种情况的，Ny 代表沿 y 方向的子晶格数目
% sublattice 的输入值有 1、2、3 三种，对应 A、B、C 三种原子
% type 的输入值有 1、2、3、4 四种，对应 AB-C \ C-AB \ AB-AB \ C-C 四种类型


freedom = 4;

if type == 1 || type == 3

    if sublattice == 1
        s = ((j-1) * (3*Nx-1) + (i-1)*2) * freedom;
    elseif sublattice == 2
        s = ((j-1) * (3*Nx-1) + (i-1)*2 + 1) * freedom;
    elseif sublattice == 3
        s = ((j-1) * (3*Nx-1) + (i-1) + (2*Nx-1)) * freedom;
    end

elseif type == 2 || type == 4

    if sublattice == 1
        s = ((j-1) * (3*Nx-1) + (i-1)*2 + Nx) * freedom;
    elseif sublattice == 2
        s = ((j-1) * (3*Nx-1) + (i-1)*2 + Nx + 1) * freedom;
    elseif sublattice == 3
        s = ((j-1) * (3*Nx-1) + (i-1)) * freedom;
    end

end



