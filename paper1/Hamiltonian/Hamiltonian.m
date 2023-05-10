function[H] = Hamiltonian(k)

t = 1;
lambda_SO = 0.1; % Fermi Surface 附近打开一个能隙
lambda_v  = 0.6; % Fermi Surface 附近上下自旋劈裂
gamma = 1;     % 粒子空穴对称性破缺
Mx = 0.0;
My = 0.0; 
M = [Mx, My, 0]; % 不同自旋取向的能带上下劈裂开

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

H = zeros(4, 4);

% kx = 0.3;
% ky = 0.4;
% k = [kx, ky];

H0 = zeros(2, 2);
H0(1, 2) = t*(exp( 1i*(dot(k,a1))) + exp( 1i*(dot(k,a2))) + exp( 1i*(dot(k,a3))));
H0(2, 1) = t*(exp(-1i*(dot(k,a1))) + exp(-1i*(dot(k,a2))) + exp(-1i*(dot(k,a3))));

H_SO = 2*lambda_SO*(sin(dot(k,b1))+sin(dot(k,b2))+sin(dot(k,b3))).*sigmaz;

H_v = lambda_v.*sigmaz;

H = [ H0 + H_SO + H_v, zeros(2,2); zeros(2,2), H0 - H_SO + H_v];

sigmax0 = [zeros(2,2),    sigma0;     sigma0, zeros(2,2)];
sigmay0 = [zeros(2,2), 1i*sigma0; -1i*sigma0, zeros(2,2)];
sigmaxz = [zeros(2,2),    sigmaz;     sigmaz, zeros(2,2)];
sigmayz = [zeros(2,2), 1i*sigmaz; -1i*sigmaz, zeros(2,2)];

H_M = (1+gamma)/2 * (Mx*sigmax0 + My*sigmay0) + (1-gamma)/2 * (Mx*sigmaxz + My*sigmayz);

H = H + H_M;


end

