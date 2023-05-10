function[H] = OBC_Hamiltonian(k,n,type)

t = 1;
lambda_SO = 0.1;
lambda_v  = 0.0;
gamma = 1;
Mx = 0.0;
My = 0.0;

if type == 2 || type == 3 || type == 4
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

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 针对六角晶格 Honeycomb，并考虑开边界为 y 方向，x 方向为周期结构
    % 开边界端原胞为 A-A 锯齿结构
    % 给定沿 y 方向有 n 个 Honeycomb 结构，该方向将包含 A-B 结构数 N = 2n

    N = n;

    H = zeros(4*N, 4*N);

    % 第一项（Hopping）

    H0 = zeros(4*N, 4*N);
    H0_mid = kron(eye(N),kron(sigma0,[0, t*(exp(1i*(dot(k,a2))) + exp(1i*(dot(k,a3)))); 0, 0]));
    H0_hopping = circshift(eye(N),-1);
    H0_hopping(N,1) = 0;
    H0_hopping = kron(H0_hopping,kron(sigma0,[0, 0; t*exp(-1i*(dot(k,a1))), 0]));
    H0 = H0_mid + H0_hopping;
    H0 = H0 + H0';

    % 第二项（SO）

    H_SO = zeros(4*N, 4*N);

    H_SO_mid = kron(eye(N),kron(sigmaz,sigmaz*(exp(1i*dot(k,b1))-exp(-1i*dot(k,b1)))));
    H_SO_up = circshift(eye(N),-1);
    H_SO_up(N,1) = 0;
    H_SO_up = kron(H_SO_up,kron(sigmaz,sigmaz*(exp(1i*dot(k,b2))-exp(-1i*dot(k,b3)))));
    H_SO_down = -H_SO_up';

    H_SO = 1i*lambda_SO*(H_SO_mid + H_SO_up + H_SO_down);

    % 第三项（v）

    H_v = zeros(4*N, 4*N);

    H_v = kron(eye(N),kron(sigma0,kron(sigmaz,lambda_v)));

    % 第四项

    M = [1, 0; 0, gamma];
    H_M0 = (Mx*kron(sigmax, M) + My*kron(sigmay, M));
    H_M = kron(eye(N), H_M0);

    % 总的 Hamiltonian

    H = H0 + H_SO + H_v + H_M;

    if type == 3
        H(4*(N-1)+4,:) = [];
        H(4*(N-1)+2,:) = [];
        H(:,4*(N-1)+4) = [];
        H(:,4*(N-1)+2) = [];
        H(3,:) = [];
        H(1,:) = [];
        H(:,3) = [];
        H(:,1) = [];
    elseif type == 4
        H(4*(N-1)+4,:) = [];
        H(4*(N-1)+2,:) = [];
        H(:,4*(N-1)+4) = [];
        H(:,4*(N-1)+2) = [];
    end
elseif type == 1
    
    a1 = [    1,          0];
    a2 = [ -1/2, -sqrt(3)/2];
    a3 = [ -1/2,  sqrt(3)/2];

    b1 = a2 - a3;
    b2 = a3 - a1;
    b3 = a1 - a2;

    sigma0 = [ 1,  0;   0,  1];
    sigmax = [ 0,  1;   1,  0];
    sigmay = [ 0, 1i; -1i,  0];
    sigmaz = [ 1,  0;   0, -1];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 针对六角晶格 Honeycomb，并考虑开边界为 y 方向，x 方向为周期结构
    % 开边界端原胞为 A-A 锯齿结构
    % 给定沿 y 方向有 n 个 Honeycomb 结构，该方向将包含 A-B 结构数 N = 2n

    N = n;

    H = zeros(4*N, 4*N);

    % 第一项（Hopping）

    H0 = zeros(4*N, 4*N);
    H0_mid  = kron(eye(N), kron(sigma0,[0, t*exp( 1i*(dot(k,a1))); 0, 0]));
    H0_up   = circshift(eye(N),-1);
    H0_up(N,1) = 0;
    H0_down = H0_up';
    H0_up   = kron(H0_up,  kron(sigma0,[0, t*exp( 1i*(dot(k,a3))); 0, 0]));
    H0_down = kron(H0_down,kron(sigma0,[0, t*exp( 1i*(dot(k,a2))); 0, 0]));
    H0 = H0_mid + H0_up + H0_down;
    H0 = H0 + H0';

    % 第二项（SO）

    H_SO = zeros(4*N, 4*N);
    H_SO_up = circshift(eye(N),-1);
    H_SO_up(N,1) = 0;
    H_SO_upup = circshift(H_SO_up,-1);
    H_SO_upup(N,2) = 0;
    
    H_SO_up = kron(H_SO_up,kron(sigmaz,sigmaz*(exp(1i*dot(k,b2))+exp(1i*dot(k,b3)))));
    H_SO_down = -H_SO_up';
    H_SO_upup = kron(H_SO_upup,kron(sigmaz,sigmaz*(-exp(-1i*dot(k,b1)))));
    H_SO_downdown = -H_SO_upup';
    
    H_SO = 1i*lambda_SO*(H_SO_up + H_SO_down + H_SO_upup + H_SO_downdown);

    % 第三项（v）

    H_v = zeros(4*N, 4*N);

    H_v = kron(eye(N),kron(sigma0,kron(sigmaz,lambda_v)));

    % 第四项

    M = [1, 0; 0, gamma];
    H_M0 = (Mx*kron(sigmax, M) + My*kron(sigmay, M));
    H_M = kron(eye(N), H_M0);

    % 总的 Hamiltonian

    H = H0 + H_SO + H_v + H_M;
end

    
end




