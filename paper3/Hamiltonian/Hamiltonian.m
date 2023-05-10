function[H] = Hamiltonian(k)

t = 1;
lambda_SO = 0.1; % Fermi Surface 附近打开一个能隙
mu  = 0.1;       % Fermi Surface 附近上下自旋劈裂
Delta0 = 0.2;    %  on-site 超导配对项
Delta1 = 0.2;    %  nearest-neighbor 超导配对项
Delta2 = 0.2;    %  next-nearest-neighbor 超导配对项

a1 = [         0,    1];
a2 = [ sqrt(3)/2, -1/2];
a3 = [-sqrt(3)/2, -1/2];
b1 = a2 - a3;
b2 = a3 - a1;
b3 = a1 - a2;

sigma0 = [ 1,  0;   0,  1];
sigmax = [ 0,  1;   1,  0];
sigmay = [ 0, 1i; -1i,  0];
sigmaz = [ 1,  0;   0, -1];

% 以下是理想晶格结构的哈密顿量

% kx = 0.3;
% ky = 0.4;
% k = [kx, ky];

H   = t*(cos(dot(k,a1))+cos(dot(k,a2))+cos(dot(k,a3)))*kron(sigmaz,kron(sigma0,sigmax)) ...
    + t*(sin(dot(k,a1))+sin(dot(k,a2))+sin(dot(k,a3)))*kron(sigmaz,kron(sigma0,sigmay)) ...
    + 2*lambda_SO*(sin(dot(k,b1))+sin(dot(k,b2))+sin(dot(k,b3)))*kron(sigma0,kron(sigmaz,sigmaz)) ...
    - Delta1*(cos(dot(k,a1))+cos(dot(k,a2))+cos(dot(k,a3)))*kron(sigmay,kron(sigmay,sigmax)) ...
    - Delta1*(sin(dot(k,a1))+sin(dot(k,a2))+sin(dot(k,a3)))*kron(sigmay,kron(sigmay,sigmay)) ...
    - Delta0*kron(sigmay,kron(sigmay,sigma0)) ...
    - Delta2*2*(cos(dot(k,b1))+cos(dot(k,b2))+cos(dot(k,b3)))*kron(sigmay,kron(sigmay,sigma0)) ...
    - mu*kron(sigmaz,kron(sigma0,sigma0));

end

