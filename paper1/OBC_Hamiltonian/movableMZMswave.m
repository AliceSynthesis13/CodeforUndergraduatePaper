%% movableMZM (x方向为周期边界条件(圆柱面),y方向为开边界)
clear all;
clc;
tot = 0;
global Nx;
global Ny;
global Nmovs;
global Nmov;
Nx = 20;
Ny = 20;
% 以下是用来在最上面一层添加sublattice的参数
Nmovs = 1; % Nmovs为起始添加sublattice的位置
Nmov  = 0; % Nmov为添加的sublattice的个数, Nmov < Ny以便观察到现象

t = 1;
ilso = 1i*0.1; % ilso=-ilso;%后面的ilambdaso项v_ij取反.
Deltad = 0.0;
Delta0 = 0.0;
mu = 0.0;
B = 0.0;
gamma = 0.3;
Mx = 0.2;
My = 0.2;

N0 = 0; % N0决定第一行是什么sublattice,0则为B,1则为A.(N0只能为0或1)

braodong=0; Nbrao=0; % beard边界加local的扰动
zraodong=0; Nzrao=0; % zigzag边界加local的扰动

%输出画图参数
r=3.5;
pi=2*3.14159;
l=3;
csize=100;
sensitive=1;
linewidth=0.5;
circlewidth=0.3;

spectrum=0; % 是否输出能谱(0为不输出,1为输出)
minex=0.26; % 画图时的最小激发能
cylinder=0; % 是否输出圆柱面的波函数分布
plane=1;    % 是否输出平面图
lattice=0;  % 是否输出晶格图(在plane=1的情况下)

%第一项最近邻相互作用
for i=1:1:Nx
    for j=1:1:Ny
        if (mod(i,4)==mod(N0+1,4))
            for k=1:1:2
                if i>1
                    matrix(posi(i,j)+k,posi(i-1,j)+k)=t;
                    matrix(posi(i-1,j)+k+2,posi(i,j)+k+2)=-t;
                end
                if i<Nx
                    matrix(posi(i,j)+k,posi(i+1,j)+k)=t;
                    matrix(posi(i+1,j)+k+2,posi(i,j)+k+2)=-t;
                    matrix(posi(i,j)+k,posi(i+1,mod1(j+1))+k)=t;
                    matrix(posi(i+1,mod1(j+1))+k+2,posi(i,j)+k+2)=-t;
                end
            end
        end
        if (mod(i,4)==mod(N0+2,4))
            for k=1:1:2
                if i>1
                    matrix(posi(i,j)+k,posi(i-1,j)+k)=t;
                    matrix(posi(i-1,j)+k+2,posi(i,j)+k+2)=-t;
                    matrix(posi(i,j)+k,posi(i-1,mod1(j-1))+k)=t;
                    matrix(posi(i-1,mod1(j-1))+k+2,posi(i,j)+k+2)=-t;
                end
                if i<Nx
                    matrix(posi(i,j)+k,posi(i+1,j)+k)=t;
                    matrix(posi(i+1,j)+k+2,posi(i,j)+k+2)=-t;
                end
            end
        end
        if (mod(i,4)==mod(N0+3,4))
            for k=1:1:2
                if i>1
                    matrix(posi(i,j)+k,posi(i-1,j)+k)=t;
                    matrix(posi(i-1,j)+k+2,posi(i,j)+k+2)=-t;
                end
                if i<Nx
                    matrix(posi(i,j)+k,posi(i+1,j)+k)=t;
                    matrix(posi(i+1,j)+k+2,posi(i,j)+k+2)=-t;
                    matrix(posi(i,j)+k,posi(i+1,mod1(j-1))+k)=t;
                    matrix(posi(i+1,mod1(j-1))+k+2,posi(i,j)+k+2)=-t;
                end
            end
        end
        if (mod(i,4)==mod(N0+4,4))
            for k=1:1:2
                if i>1
                    matrix(posi(i,j)+k,posi(i-1,j)+k)=t;
                    matrix(posi(i-1,j)+k+2,posi(i,j)+k+2)=-t;
                    matrix(posi(i,j)+k,posi(i-1,mod1(j+1))+k)=t;
                    matrix(posi(i-1,mod1(j+1))+k+2,posi(i,j)+k+2)=-t;
                end
                if i<Nx
                    matrix(posi(i,j)+k,posi(i+1,j)+k)=t;
                    matrix(posi(i+1,j)+k+2,posi(i,j)+k+2)=-t;
                end
            end
        end
    end
end

%第二项次近邻相互作用
for i=1:1:Nx
    for j=1:1:Ny
        if  (mod(i,4)==mod(N0+1,4))
            for k=1:1:2
                matrix(posi(i,j)+k,posi(i,mod1(j+1))+k)=ilso*power(-1,k-1);
                matrix(posi(i,mod1(j+1))+k+2,posi(i,j)+k+2)=-ilso*power(-1,k-1);
                matrix(posi(i,j)+k,posi(i,mod1(j-1))+k)=-ilso*power(-1,k-1);
                matrix(posi(i,mod1(j-1))+k+2,posi(i,j)+k+2)=ilso*power(-1,k-1);
                if i>2
                    matrix(posi(i,j)+k,posi(i-2,j)+k)=ilso*power(-1,k-1);
                    matrix(posi(i-2,j)+k+2,posi(i,j)+k+2)=-ilso*power(-1,k-1);
                    matrix(posi(i,j)+k,posi(i-2,mod1(j+1))+k)=-ilso*power(-1,k-1);
                    matrix(posi(i-2,mod1(j+1))+k+2,posi(i,j)+k+2)=ilso*power(-1,k-1);
                end
                if i<Nx-1
                    matrix(posi(i,j)+k,posi(i+2,j)+k)=ilso*power(-1,k-1);
                    matrix(posi(i+2,j)+k+2,posi(i,j)+k+2)=-ilso*power(-1,k-1);
                    matrix(posi(i,j)+k,posi(i+2,mod1(j+1))+k)=-ilso*power(-1,k-1);
                    matrix(posi(i+2,mod1(j+1))+k+2,posi(i,j)+k+2)=ilso*power(-1,k-1);
                end
            end
        end
        if (mod(i,4)==mod(N0+2,4))
            for k=1:1:2
                matrix(posi(i,j)+k,posi(i,mod1(j+1))+k)=-ilso*power(-1,k-1);
                matrix(posi(i,mod1(j+1))+k+2,posi(i,j)+k+2)=ilso*power(-1,k-1);
                matrix(posi(i,j)+k,posi(i,mod1(j-1))+k)=ilso*power(-1,k-1);
                matrix(posi(i,mod1(j-1))+k+2,posi(i,j)+k+2)=-ilso*power(-1,k-1);
                if i>2
                    matrix(posi(i,j)+k,posi(i-2,j)+k)=ilso*power(-1,k-1);
                    matrix(posi(i-2,j)+k+2,posi(i,j)+k+2)=-ilso*power(-1,k-1);
                    matrix(posi(i,j)+k,posi(i-2,mod1(j-1))+k)=-ilso*power(-1,k-1);
                    matrix(posi(i-2,mod1(j-1))+k+2,posi(i,j)+k+2)=ilso*power(-1,k-1);
                end
                if i<Nx-1
                    matrix(posi(i,j)+k,posi(i+2,j)+k)=ilso*power(-1,k-1);
                    matrix(posi(i+2,j)+k+2,posi(i,j)+k+2)=-ilso*power(-1,k-1);
                    matrix(posi(i,j)+k,posi(i+2,mod1(j-1))+k)=-ilso*power(-1,k-1);
                    matrix(posi(i+2,mod1(j-1))+k+2,posi(i,j)+k+2)=ilso*power(-1,k-1);
                end
            end
        end
        if (mod(i,4)==mod(N0+3,4))
            for k=1:1:2
                matrix(posi(i,j)+k,posi(i,mod1(j+1))+k)=ilso*power(-1,k-1);
                matrix(posi(i,mod1(j+1))+k+2,posi(i,j)+k+2)=-ilso*power(-1,k-1);
                matrix(posi(i,j)+k,posi(i,mod1(j-1))+k)=-ilso*power(-1,k-1);
                matrix(posi(i,mod1(j-1))+k+2,posi(i,j)+k+2)=ilso*power(-1,k-1);
                if i>2
                    matrix(posi(i,j)+k,posi(i-2,j)+k)=-ilso*power(-1,k-1);
                    matrix(posi(i-2,j)+k+2,posi(i,j)+k+2)=ilso*power(-1,k-1);
                    matrix(posi(i,j)+k,posi(i-2,mod1(j-1))+k)=ilso*power(-1,k-1);
                    matrix(posi(i-2,mod1(j-1))+k+2,posi(i,j)+k+2)=-ilso*power(-1,k-1);
                end
                if i<Nx-1
                    matrix(posi(i,j)+k,posi(i+2,j)+k)=-ilso*power(-1,k-1);
                    matrix(posi(i+2,j)+k+2,posi(i,j)+k+2)=ilso*power(-1,k-1);
                    matrix(posi(i,j)+k,posi(i+2,mod1(j-1))+k)=ilso*power(-1,k-1);
                    matrix(posi(i+2,mod1(j-1))+k+2,posi(i,j)+k+2)=-ilso*power(-1,k-1);
                end
            end
        end
        if (mod(i,4)==mod(N0+4,4))
            for k=1:1:2
                matrix(posi(i,j)+k,posi(i,mod1(j+1))+k)=-ilso*power(-1,k-1);
                matrix(posi(i,mod1(j+1))+k+2,posi(i,j)+k+2)=ilso*power(-1,k-1);
                matrix(posi(i,j)+k,posi(i,mod1(j-1))+k)=ilso*power(-1,k-1);
                matrix(posi(i,mod1(j-1))+k+2,posi(i,j)+k+2)=-ilso*power(-1,k-1);
                if i>2
                    matrix(posi(i,j)+k,posi(i-2,j)+k)=-ilso*power(-1,k-1);
                    matrix(posi(i-2,j)+k+2,posi(i,j)+k+2)=ilso*power(-1,k-1);
                    matrix(posi(i,j)+k,posi(i-2,mod1(j+1))+k)=ilso*power(-1,k-1);
                    matrix(posi(i-2,mod1(j+1))+k+2,posi(i,j)+k+2)=-ilso*power(-1,k-1);
                end
                if i<Nx-1
                    matrix(posi(i,j)+k,posi(i+2,j)+k)=-ilso*power(-1,k-1);
                    matrix(posi(i+2,j)+k+2,posi(i,j)+k+2)=ilso*power(-1,k-1);
                    matrix(posi(i,j)+k,posi(i+2,mod1(j+1))+k)=ilso*power(-1,k-1);
                    matrix(posi(i+2,mod1(j+1))+k+2,posi(i,j)+k+2)=-ilso*power(-1,k-1);
                end
            end
        end
    end
end

% s-wave超导配对
for i=1:1:Nx
    for j=1:1:Ny
        matrix(posi(i,j)+1,posi(i,j)+4)=Delta0;
        matrix(posi(i,j)+4,posi(i,j)+1)=Delta0;
        matrix(posi(i,j)+2,posi(i,j)+3)=-Delta0;
        matrix(posi(i,j)+3,posi(i,j)+2)=-Delta0;
    end
end
% 次近邻超导配对作用
for i=1:1:Nx
    for j=1:1:Ny
        matrix(posi(i,j)+1,posi(i,mod1(j+1))+4)=Deltad;
        matrix(posi(i,mod1(j+1))+4,posi(i,j)+1)=Deltad;
        matrix(posi(i,mod1(j+1))+2,posi(i,j)+3)=-Deltad;
        matrix(posi(i,j)+3,posi(i,mod1(j+1))+2)=-Deltad;
        matrix(posi(i,j)+1,posi(i,mod1(j-1))+4)=Deltad;
        matrix(posi(i,mod1(j-1))+4,posi(i,j)+1)=Deltad;
        matrix(posi(i,mod1(j-1))+2,posi(i,j)+3)=-Deltad;
        matrix(posi(i,j)+3,posi(i,mod1(j-1))+2)=-Deltad;
        if (mod(i,4)==mod(N0+1,4)||mod(i,4)==mod(N0+4,4))
            if i>2
                matrix(posi(i,j)+1,posi(i-2,j)+4)=Deltad;
                matrix(posi(i-2,j)+4,posi(i,j)+1)=Deltad;
                matrix(posi(i-2,j)+2,posi(i,j)+3)=-Deltad;
                matrix(posi(i,j)+3,posi(i-2,j)+2)=-Deltad;
                matrix(posi(i,j)+1,posi(i-2,mod1(j+1))+4)=Deltad;
                matrix(posi(i-2,mod1(j+1))+4,posi(i,j)+1)=Deltad;
                matrix(posi(i-2,mod1(j+1))+2,posi(i,j)+3)=-Deltad;
                matrix(posi(i,j)+3,posi(i-2,mod1(j+1))+2)=-Deltad;
            end
            if i<Nx-1
                matrix(posi(i,j)+1,posi(i+2,j)+4)=Deltad;
                matrix(posi(i+2,j)+4,posi(i,j)+1)=Deltad;
                matrix(posi(i+2,j)+2,posi(i,j)+3)=-Deltad;
                matrix(posi(i,j)+3,posi(i+2,j)+2)=-Deltad;
                matrix(posi(i,j)+1,posi(i+2,mod1(j+1))+4)=Deltad;
                matrix(posi(i+2,mod1(j+1))+4,posi(i,j)+1)=Deltad;
                matrix(posi(i+2,mod1(j+1))+2,posi(i,j)+3)=-Deltad;
                matrix(posi(i,j)+3,posi(i+2,mod1(j+1))+2)=-Deltad;
            end
        end
        if (mod(i,4)==mod(N0+2,4)||mod(i,4)==mod(N0+3,4))
            if i>2
                matrix(posi(i,j)+1,posi(i-2,j)+4)=Deltad;
                matrix(posi(i-2,j)+4,posi(i,j)+1)=Deltad;
                matrix(posi(i-2,j)+2,posi(i,j)+3)=-Deltad;
                matrix(posi(i,j)+3,posi(i-2,j)+2)=-Deltad;
                matrix(posi(i,j)+1,posi(i-2,mod1(j-1))+4)=Deltad;
                matrix(posi(i-2,mod1(j-1))+4,posi(i,j)+1)=Deltad;
                matrix(posi(i-2,mod1(j-1))+2,posi(i,j)+3)=-Deltad;
                matrix(posi(i,j)+3,posi(i-2,mod1(j-1))+2)=-Deltad;
            end
            if i<Nx-1
                matrix(posi(i,j)+1,posi(i+2,j)+4)=Deltad;
                matrix(posi(i+2,j)+4,posi(i,j)+1)=Deltad;
                matrix(posi(i+2,j)+2,posi(i,j)+3)=-Deltad;
                matrix(posi(i,j)+3,posi(i+2,j)+2)=-Deltad;
                matrix(posi(i,j)+1,posi(i+2,mod1(j-1))+4)=Deltad;
                matrix(posi(i+2,mod1(j-1))+4,posi(i,j)+1)=Deltad;
                matrix(posi(i+2,mod1(j-1))+2,posi(i,j)+3)=-Deltad;
                matrix(posi(i,j)+3,posi(i+2,mod1(j-1))+2)=-Deltad;
            end
        end
    end
end

% 化学势
for i=1:1:Nx
    for j=1:1:Ny
        for k=1:1:2
            matrix(posi(i,j)+k,posi(i,j)+k)=matrix(posi(i,j)+k,posi(i,j)+k)-mu;
            matrix(posi(i,j)+k+2,posi(i,j)+k+2)=matrix(posi(i,j)+k+2,posi(i,j)+k+2)+mu;
        end
    end
end
% 塞曼场
for i=1:1:Nx
    for j=1:1:Ny
        matrix(posi(i,j)+1,posi(i,j)+2)=matrix(posi(i,j)+1,posi(i,j)+2)+B;
        matrix(posi(i,j)+2,posi(i,j)+1)=matrix(posi(i,j)+2,posi(i,j)+1)+B;
        matrix(posi(i,j)+3,posi(i,j)+4)=matrix(posi(i,j)+3,posi(i,j)+4)-B;
        matrix(posi(i,j)+4,posi(i,j)+3)=matrix(posi(i,j)+4,posi(i,j)+3)-B;
    end
end

% 磁场
for i=1:1:Nx
    for j=1:1:Ny
        if N0 == 0
            if mod(i,2) == 1 % A类 
                matrix(posi(i,j)+1,posi(i,j)+2) = matrix(posi(i,j)+1,posi(i,j)+2) + Mx + 1i*My;
                matrix(posi(i,j)+2,posi(i,j)+1) = matrix(posi(i,j)+2,posi(i,j)+1) + Mx - 1i*My;
                matrix(posi(i,j)+3,posi(i,j)+4) = matrix(posi(i,j)+3,posi(i,j)+4) -(Mx - 1i*My);
                matrix(posi(i,j)+4,posi(i,j)+3) = matrix(posi(i,j)+4,posi(i,j)+3) -(Mx + 1i*My);
            elseif mod(i,2) == 0 % B类
                matrix(posi(i,j)+1,posi(i,j)+2) = matrix(posi(i,j)+1,posi(i,j)+2) + gamma*(Mx + 1i*My);
                matrix(posi(i,j)+2,posi(i,j)+1) = matrix(posi(i,j)+2,posi(i,j)+1) + gamma*(Mx - 1i*My);
                matrix(posi(i,j)+3,posi(i,j)+4) = matrix(posi(i,j)+3,posi(i,j)+4) - gamma*(Mx - 1i*My);
                matrix(posi(i,j)+4,posi(i,j)+3) = matrix(posi(i,j)+4,posi(i,j)+3) - gamma*(Mx + 1i*My);
            end
        elseif N0 == 1
            if mod(i,2) == 0 % A类 
                matrix(posi(i,j)+1,posi(i,j)+2) = matrix(posi(i,j)+1,posi(i,j)+2) + Mx + 1i*My;
                matrix(posi(i,j)+2,posi(i,j)+1) = matrix(posi(i,j)+2,posi(i,j)+1) + Mx - 1i*My;
                matrix(posi(i,j)+3,posi(i,j)+4) = matrix(posi(i,j)+3,posi(i,j)+4) -(Mx - 1i*My);
                matrix(posi(i,j)+4,posi(i,j)+3) = matrix(posi(i,j)+4,posi(i,j)+3) -(Mx + 1i*My);
            elseif mod(i,2) == 1 % B类
                matrix(posi(i,j)+1,posi(i,j)+2) = matrix(posi(i,j)+1,posi(i,j)+2) + gamma*(Mx + 1i*My);
                matrix(posi(i,j)+2,posi(i,j)+1) = matrix(posi(i,j)+2,posi(i,j)+1) + gamma*(Mx - 1i*My);
                matrix(posi(i,j)+3,posi(i,j)+4) = matrix(posi(i,j)+3,posi(i,j)+4) - gamma*(Mx - 1i*My);
                matrix(posi(i,j)+4,posi(i,j)+3) = matrix(posi(i,j)+4,posi(i,j)+3) - gamma*(Mx + 1i*My);
            end
        end
    end
end

%% 在六角晶格的上方添加sublattice
if (Nmovs>=1&&Nmovs<=Ny&&Nmov<=Ny&&Nmov>0)
    %最近邻相互作用
    for j=Nmovs:1:Nmov+Nmovs-1
        if (N0==0)
            for k=1:1:2
                matrix(posi(1,j)+k,posim(j)+k)=t;
                matrix(posim(j)+k,posi(1,j)+k)=t;
                matrix(posim(j)+k+2,posi(1,j)+k+2)=-t;
                matrix(posi(1,j)+k+2,posim(j)+k+2)=-t;
            end
        end
        if (N0==1)
            for k=1:1:2
                matrix(posi(1,j)+k,posim(j)+k)=t;
                matrix(posim(j)+k,posi(1,j)+k)=t;
                matrix(posim(j)+k+2,posi(1,j)+k+2)=-t;
                matrix(posi(1,j)+k+2,posim(j)+k+2)=-t;
                matrix(posim(j)+k,posi(1,mod1(j+1))+k)=t;
                matrix(posi(1,mod1(j+1))+k,posim(j)+k)=t;
                matrix(posi(1,mod1(j+1))+k+2,posim(j)+k+2)=-t;
                matrix(posim(j)+k+2,posi(1,mod1(j+1))+k+2)=-t;
            end
        end
    end

    %次近邻相互作用
    for j=Nmovs:1:Nmov+Nmovs-1
        j=mod1(j);
        if (N0==0)
            for k=1:1:2
                if Nmov~=Ny
                    if j~=Nmovs
                        matrix(posim(j)+k,posim(mod1(j-1))+k)=ilso*power(-1,k-1);
                        matrix(posim(mod1(j-1))+k+2,posim(j)+k+2)=-ilso*power(-1,k-1);
                    end
                    if j~=mod1(Nmov+Nmovs-1)
                        matrix(posim(j)+k,posim(mod1(j+1))+k)=-ilso*power(-1,k-1);
                        matrix(posim(mod1(j+1))+k+2,posim(j)+k+2)=ilso*power(-1,k-1);
                    end
                else
                    matrix(posim(j)+k,posim(mod1(j-1))+k)=ilso*power(-1,k-1);
                    matrix(posim(mod1(j-1))+k+2,posim(j)+k+2)=-ilso*power(-1,k-1);
                    matrix(posim(j)+k,posim(mod1(j+1))+k)=-ilso*power(-1,k-1);
                    matrix(posim(mod1(j+1))+k+2,posim(j)+k+2)=ilso*power(-1,k-1);
                end
                matrix(posim(j)+k,posi(2,j)+k)=-ilso*power(-1,k-1);
                matrix(posi(2,j)+k,posim(j)+k)=ilso*power(-1,k-1);
                matrix(posi(2,j)+k+2,posim(j)+k+2)=ilso*power(-1,k-1);
                matrix(posim(j)+k+2,posi(2,j)+k+2)=-ilso*power(-1,k-1);
                matrix(posim(j)+k,posi(2,mod1(j+1))+k)=ilso*power(-1,k-1);
                matrix(posi(2,mod1(j+1))+k,posim(j)+k)=-ilso*power(-1,k-1);
                matrix(posi(2,mod1(j+1))+k+2,posim(j)+k+2)=-ilso*power(-1,k-1);
                matrix(posim(j)+k+2,posi(2,mod1(j+1))+k+2)=ilso*power(-1,k-1);
            end
        end
        if (N0==1)
            for k=1:1:2
                if j~=Nmovs
                    matrix(posim(j)+k,posim(mod1(j-1))+k)=-ilso*power(-1,k-1);
                    matrix(posim(mod1(j-1))+k+2,posim(j)+k+2)=ilso*power(-1,k-1);
                end
                if j~=mod1(Nmov+Nmovs-1)
                    matrix(posim(j)+k,posim(mod1(j+1))+k)=ilso*power(-1,k-1);
                    matrix(posim(mod1(j+1))+k+2,posim(j)+k+2)=-ilso*power(-1,k-1);
                end
                matrix(posim(j)+k,posi(2,j)+k)=ilso*power(-1,k-1);
                matrix(posi(2,j)+k,posim(j)+k)=-ilso*power(-1,k-1);
                matrix(posi(2,j)+k+2,posim(j)+k+2)=-ilso*power(-1,k-1);
                matrix(posim(j)+k+2,posi(2,j)+k+2)=ilso*power(-1,k-1);
                matrix(posim(j)+k,posi(2,mod1(j+1))+k)=-ilso*power(-1,k-1);
                matrix(posi(2,mod1(j+1))+k,posim(j)+k)=ilso*power(-1,k-1);
                matrix(posi(2,mod1(j+1))+k+2,posim(j)+k+2)=ilso*power(-1,k-1);
                matrix(posim(j)+k+2,posi(2,mod1(j+1))+k+2)=-ilso*power(-1,k-1);
            end
        end
    end

    %s-wave超导
    for j=Nmovs:1:Nmovs+Nmov-1
        matrix(posim(j)+1,posim(j)+4)=Delta0;
        matrix(posim(j)+4,posim(j)+1)=Delta0;
        matrix(posim(j)+2,posim(j)+3)=-Delta0;
        matrix(posim(j)+3,posim(j)+2)=-Delta0;
    end

    %超导次近邻相互作用
    for j=Nmovs:1:Nmov+Nmovs-1
        j=mod1(j);
        if (N0==0||N0==1)
            if Nmov~=Ny
                if j~=Nmovs
                    matrix(posim(j)+1,posim(mod1(j-1))+4)=Deltad;
                    matrix(posim(mod1(j-1))+4,posim(j)+1)=Deltad;
                    matrix(posim(mod1(j-1))+2,posim(j)+3)=-Deltad;
                    matrix(posim(j)+3,posim(mod1(j-1))+2)=-Deltad;
                end
                if j~=mod1(Nmov+Nmovs-1)
                    matrix(posim(j)+1,posim(mod1(j+1))+4)=Deltad;
                    matrix(posim(mod1(j+1))+4,posim(j)+1)=Deltad;
                    matrix(posim(mod1(j+1))+2,posim(j)+3)=-Deltad;
                    matrix(posim(j)+3,posim(mod1(j+1))+2)=-Deltad;
                end
            else
                matrix(posim(j)+1,posim(mod1(j-1))+4)=Deltad;
                matrix(posim(mod1(j-1))+4,posim(j)+1)=Deltad;
                matrix(posim(mod1(j-1))+2,posim(j)+3)=-Deltad;
                matrix(posim(j)+3,posim(mod1(j-1))+2)=-Deltad;
                matrix(posim(j)+1,posim(mod1(j+1))+4)=Deltad;
                matrix(posim(mod1(j+1))+4,posim(j)+1)=Deltad;
                matrix(posim(mod1(j+1))+2,posim(j)+3)=-Deltad;
                matrix(posim(j)+3,posim(mod1(j+1))+2)=-Deltad;
            end
            matrix(posim(j)+1,posi(2,j)+4)=Deltad;
            matrix(posi(2,j)+4,posim(j)+1)=Deltad;
            matrix(posi(2,j)+2,posim(j)+3)=-Deltad;
            matrix(posim(j)+3,posi(2,j)+2)=-Deltad;
            matrix(posi(2,j)+1,posim(j)+4)=Deltad;
            matrix(posim(j)+4,posi(2,j)+1)=Deltad;
            matrix(posim(j)+2,posi(2,j)+3)=-Deltad;
            matrix(posi(2,j)+3,posim(j)+2)=-Deltad;

            matrix(posim(j)+1,posi(2,mod1(j+1))+4)=Deltad;
            matrix(posi(2,mod1(j+1))+4,posim(j)+1)=Deltad;
            matrix(posi(2,mod1(j+1))+2,posim(j)+3)=-Deltad;
            matrix(posim(j)+3,posi(2,mod1(j+1))+2)=-Deltad;
            matrix(posi(2,mod1(j+1))+1,posim(j)+4)=Deltad;
            matrix(posim(j)+4,posi(2,mod1(j+1))+1)=Deltad;
            matrix(posim(j)+2,posi(2,mod1(j+1))+3)=-Deltad;
            matrix(posi(2,mod1(j+1))+3,posim(j)+2)=-Deltad;
        end
    end
    %beard local扰动项
    if (Nbrao<Nmov&& Nbrao>0)
        for i=Nmovs:1:Nmovs+Nbrao-1
            i=mod1(i);
            matrix(posim(i)+1,posim(i)+1)=matrix(posim(i)+1,posim(i)+1)+braodong;
            matrix(posim(i)+2,posim(i)+2)=matrix(posim(i)+2,posim(i)+2)+braodong;
            matrix(posim(i)+3,posim(i)+3)=matrix(posim(i)+3,posim(i)+3)-braodong;
            matrix(posim(i)+4,posim(i)+4)=matrix(posim(i)+4,posim(i)+4)-braodong;
        end
    end

    %zigzag local 扰动项
    if (Nzrao<Nmov && Nzrao>0)
        for i=Nmovs:-1:Nmovs-Nzrao+1
            i=mod1(i);
            matrix(posi(1,i)+1,posi(1,i)+1)=matrix(posi(1,i)+1,posi(1,i)+1)+zraodong;
            matrix(posi(1,i)+2,posi(1,i)+2)=matrix(posi(1,i)+2,posi(1,i)+2)+zraodong;
            matrix(posi(1,i)+3,posi(1,i)+3)=matrix(posi(1,i)+3,posi(1,i)+3)-zraodong;
            matrix(posi(1,i)+4,posi(1,i)+4)=matrix(posi(1,i)+4,posi(1,i)+4)-zraodong;
        end
    end

    %化学势
    for i=Nmovs:1:Nmov+Nmovs-1
        for j=1:1:2
            matrix(posim(i)+j,posim(i)+j)=matrix(posim(i)+j,posim(i)+j)-mu;
            matrix(posim(i)+j+2,posim(i)+j+2)=matrix(posim(i)+j+2,posim(i)+j+2)+mu;
        end
    end

    %塞曼场
    for i=Nmovs:1:Nmov+Nmovs-1
        matrix(posim(i)+1,posim(i)+2)=matrix(posim(i)+1,posim(i)+2)+B;
        matrix(posim(i)+2,posim(i)+1)=matrix(posim(i)+2,posim(i)+1)+B;
        matrix(posim(i)+3,posim(i)+4)=matrix(posim(i)+3,posim(i)+4)-B;
        matrix(posim(i)+4,posim(i)+3)=matrix(posim(i)+4,posim(i)+3)-B;
    end
    
    % 磁场
    for j=Nmovs:1:Nmovs+Nmov-1
        if N0 == 0
            if mod(Ny,2) == 0 % A类 
                matrix(posim(j)+1,posim(j)+2) = matrix(posim(j)+1,posim(j)+2) + Mx + 1i*My;
                matrix(posim(j)+2,posim(j)+1) = matrix(posim(j)+2,posim(j)+1) + Mx - 1i*My;
                matrix(posim(j)+3,posim(j)+4) = matrix(posim(j)+3,posim(j)+4) - (Mx - 1i*My);
                matrix(posim(j)+4,posim(j)+3) = matrix(posim(j)+4,posim(j)+3) - (Mx + 1i*My);
            elseif mod(Ny,2) == 1 % B类
                matrix(posim(j)+1,posim(j)+2) = matrix(posim(j)+1,posim(j)+2) + gamma*(Mx + 1i*My);
                matrix(posim(j)+2,posim(j)+1) = matrix(posim(j)+2,posim(j)+1) + gamma*(Mx - 1i*My);
                matrix(posim(j)+3,posim(j)+4) = matrix(posim(j)+3,posim(j)+4) - gamma*(Mx - 1i*My);
                matrix(posim(j)+4,posim(j)+3) = matrix(posim(j)+4,posim(j)+3) - gamma*(Mx + 1i*My);
            end
        elseif N0 == 1
            if mod(i,2) == 1 % A类 
                matrix(posim(j)+1,posim(j)+2) = matrix(posim(j)+1,posim(j)+2) + Mx + 1i*My;
                matrix(posim(j)+2,posim(j)+1) = matrix(posim(j)+2,posim(j)+1) + Mx - 1i*My;
                matrix(posim(j)+3,posim(j)+4) = matrix(posim(j)+3,posim(j)+4) - (Mx - 1i*My);
                matrix(posim(j)+4,posim(j)+3) = matrix(posim(j)+4,posim(j)+3) - (Mx + 1i*My);
            elseif mod(i,2) == 0 % B类
                matrix(posim(j)+1,posim(j)+2) = matrix(posim(j)+1,posim(j)+2) + gamma*(Mx + 1i*My);
                matrix(posim(j)+2,posim(j)+1) = matrix(posim(j)+2,posim(j)+1) + gamma*(Mx - 1i*My);
                matrix(posim(j)+3,posim(j)+4) = matrix(posim(j)+3,posim(j)+4) - gamma*(Mx - 1i*My);
                matrix(posim(j)+4,posim(j)+3) = matrix(posim(j)+4,posim(j)+3) - gamma*(Mx + 1i*My);
            end
        end
    end
end
%% 输出
[vet,value]=eig(matrix);
[Y,I] = sort(diag(value),'descend');        % 降序，Y 是已经排序好的数组，而 I 是 Y 每个数在原来数组（矩阵）中的位置
min(abs(Y))
Y(length(Y)/2-3:length(Y)/2+4)

if cylinder==1 % 是否输出圆柱面的波函数分布
    z0=0;dtheta=pi/Ny;
    Max=0;
    figure(1);
    for i=1:1:Nx
        if (mod(N0+1,4)==mod(i,4)||mod(N0+3,4)==mod(i,4))
            z0=z0-l;
        end
        if (mod(N0+2,4)==mod(i,4)||mod(N0+4,4)==mod(i,4))
            z0=z0-l/2;
        end
        for j=1:1:Ny
            if (mod(N0+1,4)==mod(i,4)||mod(N0+4,4)==mod(i,4))
                x(Ny*(i-1)+j)=r*cos(dtheta*(j-1));
                y(Ny*(i-1)+j)=r*sin(dtheta*(j-1));
                z(Ny*(i-1)+j)=z0;
            end
            if (mod(N0+2,4)==mod(i,4)||mod(N0+3,4)==mod(i,4))
                x(Ny*(i-1)+j)=r*cos(dtheta*(j-1)-dtheta/2);
                y(Ny*(i-1)+j)=r*sin(dtheta*(j-1)-dtheta/2);
                z(Ny*(i-1)+j)=z0;
            end
            C(Ny*(i-1)+j,2)=sensitive*abs(vet(posi(i,j)+1,2*(Nx*Ny+Nmov)));
            C(Ny*(i-1)+j,3)=C(Ny*(i-1)+j,2);
            if  C(Ny*(i-1)+j,2)>Max
                Max=C(Ny*(i-1)+j,2);
            end
            % 画晶格形状:
            if (mod(N0+1,4)==mod(i,4)||mod(N0+3,4)==mod(i,4))
                if i>1
                    plot3([x(Ny*(i-1)+j),x(Ny*(i-1)+j)],[y(Ny*(i-1)+j),y(Ny*(i-1)+j)],[z0,z0+l]);
                    hold on;
                end
            elseif (mod(N0+2,4)==mod(i,4))
                plot3([x(Ny*(i-1)+j),r*cos(dtheta*(j-1))],[y(Ny*(i-1)+j),r*sin(dtheta*(j-1))],[z0,z0+l/2],'k');
                hold on;
                plot3([x(Ny*(i-1)+j),r*cos(dtheta*(j-2))],[y(Ny*(i-1)+j),r*sin(dtheta*(j-2))],[z0,z0+l/2],'k');
                hold on;
            elseif mod(N0+4,4)==mod(i,4)
                plot3([x(Ny*(i-1)+j),r*cos(dtheta*(j-1/2))],[y(Ny*(i-1)+j),r*sin(dtheta*(j-1/2))],[z0,z0+l/2],'k');
                hold on;
                plot3([x(Ny*(i-1)+j),r*cos(dtheta*(j-3/2))],[y(Ny*(i-1)+j),r*sin(dtheta*(j-3/2))],[z0,z0+l/2],'k');
                hold on;
            end
        end
    end
    for j=1:1:Nmov
        z0=0;
        if N0==0
            z(Nx*Ny+j)=z0;
            x(Nx*Ny+j)=r*cos(dtheta*(j-2+Nmovs));
            y(Nx*Ny+j)=r*sin(dtheta*(j-2+Nmovs));
            plot3([x(Nx*Ny+j),x(Nx*Ny+j)],[y(Nx*Ny+j),y(Nx*Ny+j)],[z0,z0-l],'k');
            hold on;
        end
        if N0==1
            z(Nx*Ny+j)=z0;
            x(Nx*Ny+j)=r*cos(dtheta*(j-2+Nmovs)-dtheta/2);
            y(Nx*Ny+j)=r*sin(dtheta*(j-2+Nmovs)-dtheta/2);
            plot3([x(Ny*(i-1)+j),r*cos(dtheta*(j-2+Nmovs))],[y(Ny*(i-1)+j),r*sin(dtheta*(j-2+Nmovs))],[z0,z0-l/2],'k');
            hold on;
            plot3([x(Ny*(i-1)+j),r*cos(dtheta*(j-2+Nmovs)-dtheta)],[y(Ny*(i-1)+j),r*sin(dtheta*(j-2+Nmovs)-dtheta)],[z0,z0-l/2],'k');
            hold on;
        end
        C(Nx*Ny+j,2)=sensitive*abs(vet(4*(Nx*Ny+j-1)+1,2*(Nx*Ny+Nmov)));
        C(Nx*Ny+j,3)=C(Nx*Ny+j,2);
        if  C(Nx*Ny+j,2)>Max
            Max=C(Nx*Ny+j,2);
        end
    end
%     C(:,2)=1/Max*(Max-C(:,2));
%     C(:,3)=1/Max*(Max-C(:,3));
%     C(:,1)=1;%红色分布
%     C(:,2)=1;
%     C(:,3)=1/Max*(Max-C(:,3));
%     C(:,1)=C(:,3);%绿色分布
    C(:,1)=0;
    C(:,2)=0;C(:,3)=1;%纯画图
    scatter3(x,y,z,csize,C,'filled');
    xlabel('x');
    ylabel('y');
    zlabel('z');
end
if spectrum==1
    x=1:1:4*(Nx*Ny+Nmov);
    eigen=eig(matrix);
    figure(2);
    scatter(x,eigen);
    hold on;
    %    x=1:1:4;
    %    c1(x,1)=85/255;
    %    c1(x,2)=26/255;
    %    c1(x,3)=139/255;%紫色矩阵
    x=2*(Nx*Ny+Nmov)-1:1:2*(Nx*Ny+Nmov)+2;
    scatter(x,eigen(x,1),'r');%标记零模(无塞曼场)
    %    x=2*(Nx*Ny+Nmov):1:2*(Nx*Ny+Nmov)+1;
    %    scatter(x,eigen(x,1),'r');%标记零模(有塞曼场)
    %    scatter(x,eigen(x,1),[],c1);%标记零模(紫色)
    yticks([-3 -minex 0 minex 3])
    xticklabels({''})
    yticklabels({''})
end
if plane==1
    y0=0;
    Max=0;
    figure(3);
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
            tot=tot+bohanshu;
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
        tot=tot+bohanshu;
        C(Nx*Ny+j,2)=sensitive*bohanshu;
        C(Nx*Ny+j,3)=C(Nx*Ny+j,2);
        if  C(Nx*Ny+j,2)>Max
            Max=C(Nx*Ny+j,2);
        end
    end
    C(:,2)=1/Max*(Max-C(:,2));
    C(:,3)=1/Max*(Max-C(:,3));
    C(:,1)=1;%红色分布
    % C(:,1)=C(:,3);
    % C(:,2)=1-(255-26)/255*(C(:,2)/Max);
    % C(:,3)=1-(255-139)/255*(C(:,3)/Max);
    % C(:,1)=1-(255-85)/255*(C(:,1)/Max);%紫色分布
    % C(:,2)=1;
    % C(:,1)=1/Max*(Max-C(:,3));
    % C(:,3)=1;%Cyan色分布
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
    set(gcf,'unit','normalized','position',[0.1,0.1,0.8,0.8]);
    saveas(gcf,'.\Figure\Fig.jpg')
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