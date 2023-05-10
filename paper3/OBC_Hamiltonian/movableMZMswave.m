%% movableMZM (x方向为周期边界条件(圆柱面),y方向为开边界)
clear all;
clc;
tot = 0;
global Nx;
global Ny;
global Nmovs;
global Nmov;
Nx = 11;
Ny = 5;
% 以下是用来在最上面一层添加sublattice的参数
Nmovs = 1; % Nmovs为起始添加sublattice的位置
Nmov  = 0; % Nmov为添加的sublattice的个数, Nmov < Ny以便观察到现象

N0 = 1; % N0决定第一行是什么sublattice,0则为B,1则为A.(N0只能为0或1)

spectrum = 0; % 是否输出能谱(0为不输出,1为输出)
minex = 0.26; % 画图时的最小激发能
plane = 1;    % 是否输出平面图
lattice = 1;  % 是否输出晶格图(在plane=1的情况下)

t = 1;
ilso = 1i*0.1; % ilso=-ilso;%后面的ilambdaso项v_ij取反.
Deltad = 0.2;
Delta0 = 0.2;
mu = 0.0;
B  = 0.0;

str1 = strcat('.\Figure\OBC_M_E_type1','.jpg');
str2 = strcat('.\Figure\R_OBC_M_DoS_type3','.jpg');

braodong = 2; Nbrao = 0; % beard边界加local的扰动
zraodong = 0; Nzrao = 0; % zigzag边界加local的扰动

%输出画图参数
r=3.5;
pi=2*3.14159;
l=3;
csize=100;
sensitive=1;
linewidth=0.5;
circlewidth=0.3;


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

end
%% 输出
[vet,value]=eig(matrix);
[Y,I] = sort(diag(value),'descend');        % 降序，Y 是已经排序好的数组，而 I 是 Y 每个数在原来数组（矩阵）中的位置
minex = min(abs(Y));

if spectrum==1
    figure(2);
    scatter(0.5:1:19.5,-Y(length(Y)/2-9:length(Y)/2+10),'linewidth',1);hold on;
    scatter(8.5:1:11.5,-Y(length(Y)/2-1:length(Y)/2+2),'red','linewidth',2);
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
    saveas(gcf,str1);
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
    set(gcf,'unit','normalized','position',[0.1,0.1,0.25,0.25]);
    saveas(gcf,str2);
end


%% 自定义函数
function P=posi(i,j)
global Nx; global Ny;
P=4*(Ny*(i-1)+(j-1));
end
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