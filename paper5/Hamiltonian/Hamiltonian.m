function[H] = Hamiltonian(k)

t = 1;
lambda_SO = 0.2; % Fermi Surface 附近打开一个能隙
mu = 0.0;        % Fermi Surface 附近上下自旋劈裂
Delta = 0.0;     % 超导配对项

kx = k(1);
ky = k(2);

sigma0 = [ 1,  0;   0,  1];
sigmax = [ 0,  1;   1,  0];
sigmay = [ 0, 1i; -1i,  0];
sigmaz = [ 1,  0;   0, -1];

% 以下是理想晶格结构的哈密顿量

H = zeros(12, 12);

H0 = -2*t*[0,cos(kx),cos(ky);cos(kx),0,0;cos(ky),0,0];
H0 = kron(sigmaz,kron(sigma0,H0));

H_SO = 4*lambda_SO*[0,0,0;0,0,-1i*sin(kx)*sin(ky);0,1i*sin(kx)*sin(ky),0];
H_SO = kron(sigma0,kron(sigmaz,H_SO));

H_mu = mu*kron(sigmaz,kron(sigma0,eye(3)));

H_SC = zeros(3,3);
H_SC(1,2) =  cos(kx);
H_SC(1,3) = -cos(ky);
H_SC(2,1) =  cos(kx);
H_SC(3,1) = -cos(ky);
H_SC = Delta*kron(sigmay,kron(sigmay,H_SC));

H = H0 + H_SO - H_mu + H_SC;


end

