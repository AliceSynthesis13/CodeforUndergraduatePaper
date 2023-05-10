function[H] = OBC_Hamiltonian(k,n,type)

t = 1;
lambda_SO = 0.2;
mu = 0.8;
Delta = 0.3;

sigma0 = [ 1,  0;   0,  1];
sigmax = [ 0,  1;   1,  0];
sigmay = [ 0, 1i; -1i,  0];
sigmaz = [ 1,  0;   0, -1];

UP   = [1, 0; 0, 0];
DOWN = [0, 0; 0, 1];

N = n;

kx = k(1);
ky = k(2);


H = zeros(3*2*N*2,3*2*N*2);

% 以下是基础哈密顿量
H0 = zeros(3*2*N*2,3*2*N*2);
% 初基原胞自由度 × 上下自旋 × 原胞数 × 超导导致的电子空穴 double

% 电子部分
% A to B + h.c.
U = [0, 2*cos(kx), 0; 2*cos(kx), 0, 0; 0, 0, 0];
H0_a = kron(eye(N),kron(sigmaz,kron(sigma0,U)));
% A to down C + h.c.
V = [0, 0, 1; 0, 0, 0; 0, 0, 0];
H_up = circshift(eye(N),-1);
H_up(N,1) = 0;
H_down = H_up';
H0_b = kron(H_down,kron(sigmaz,kron(sigma0,V)));
H0_b = H0_b + H0_b';
% A to up C + h.c.
W = [0, 0, 1; 0, 0, 0; 0, 0, 0];
H0_c = kron(eye(N),kron(sigmaz,kron(sigma0,W)));
H0_c = H0_c + H0_c';
H0 = -t*(H0_a + H0_b + H0_c);

% 以下是 SOC 哈密顿量
H_SO = zeros(3*2*N*2,3*2*N*2);

% 电子部分
W = [0, 0, 0; 0, 0, 0; 0, 1, 0];
H_up = circshift(eye(N),-1);
H_up(N,1) = 0;
H_SO_a = -1i*lambda_SO*(exp(1i*kx)-exp(-1i*kx)).*kron(H_up,  kron(sigma0,kron(sigmaz,W)));
H_SO_a = H_SO_a + H_SO_a';
M = [0, 0, 0; 0, 0, 0; 0, 1, 0];
H_SO_b =  1i*lambda_SO*(exp(1i*kx)-exp(-1i*kx)).*kron(eye(N),kron(sigma0,kron(sigmaz,M)));
H_SO_b = H_SO_b + H_SO_b';
H_SO = H_SO_a + H_SO_b;

% 化学势部分
H_mu = mu*kron(eye(N),kron(sigmaz,eye(6)));

H = H0 + H_SO - H_mu;
H_normal = H;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 接下来将考虑超导态下的结果
H_up = circshift(eye(N),-1);
H_up(N,1) = 0;
H_down = H_up';

H_SC = zeros(3*2*N*2,3*2*N*2);

% 对于初基原胞的 A 原子，相互作用项为：
H_SC_AB = zeros(12,12);
H_SC_AB( 1,11) = Delta*cos(kx)*2;
H_SC_AB( 4, 8) = Delta*cos(kx)*2;
H_SC_AB = kron(eye(N),H_SC_AB);
H_SC_AB = H_SC_AB + H_SC_AB';

H_SC_AC1 = zeros(12,12);
H_SC_AC1( 1,12) = -Delta;
H_SC_AC1( 4, 9) = -Delta;
H_SC_AC1 = kron(H_down,H_SC_AC1);
H_SC_AC1 = H_SC_AC1 + H_SC_AC1';
H_SC_AC2 = zeros(12,12);
H_SC_AC2( 1,12) = -Delta;
H_SC_AC2( 4, 9) = -Delta;
H_SC_AC2 = kron(eye(N),H_SC_AC2);
H_SC_AC2 = H_SC_AC2 + H_SC_AC2';
H_SC_AC = H_SC_AC1 + H_SC_AC2;

H_SC_A = H_SC_AB + H_SC_AC;

% 对于初基原胞的 B 原子，相互作用项为：
H_SC_BA = zeros(12,12);
H_SC_BA( 2,10) = Delta*cos(kx)*2;
H_SC_BA( 5, 7) = Delta*cos(kx)*2;
H_SC_BA = kron(eye(N),H_SC_BA);
H_SC_BA = H_SC_BA + H_SC_BA';
H_SC_B = H_SC_BA;

% 对于初基原胞的 C 原子，相互作用项为：
H_SC_CA1 = zeros(12,12);
H_SC_CA1( 3,10) = -Delta;
H_SC_CA1( 6, 7) = -Delta;
H_SC_CA1 = kron(H_up,H_SC_CA1);
H_SC_CA1 = H_SC_CA1 + H_SC_CA1';
H_SC_CA2 = zeros(12,12);
H_SC_CA2( 3,10) = -Delta;
H_SC_CA2( 6, 7) = -Delta;
H_SC_CA2 = kron(eye(N),H_SC_CA2);
H_SC_CA2 = H_SC_CA2 + H_SC_CA2';
H_SC_C = H_SC_CA1 + H_SC_CA2;

% 总的超导哈密顿量为：
H_SC = H_SC_A + H_SC_B + H_SC_C;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 总的哈密顿量为：
H = H_normal + H_SC;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 当采用第二种边界条件时，解开下面

if type == 3
    H(12*(N-1)+12,:) = [];
    H(12*(N-1)+9, :) = [];
    H(12*(N-1)+6, :) = [];
    H(12*(N-1)+3, :) = [];

    H(:,12*(N-1)+12) = [];
    H(:,12*(N-1)+9 ) = [];
    H(:,12*(N-1)+6 ) = [];
    H(:,12*(N-1)+3 ) = [];
elseif type == 4
    H(11, :) = [];
    H(10, :) = [];
    H( 8, :) = [];
    H( 7, :) = [];
    H( 5, :) = [];
    H( 4, :) = [];
    H( 2, :) = [];
    H( 1, :) = [];

    H(:, 11) = [];
    H(:, 10) = [];
    H(:, 8 ) = [];
    H(:, 7 ) = [];
    H(:, 5 ) = [];
    H(:, 4 ) = [];
    H(:, 2 ) = [];
    H(:, 1 ) = [];
end

end
