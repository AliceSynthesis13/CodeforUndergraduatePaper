function[H] = OBC_Hamiltonian(k,n,type)

t = 1;
lambda_SO = 0.1;
mu  = 0.0;
Delta0 = 0.2;
Delta1 = 0.0;
Delta2 = 0.1;

Delta0 = Delta0/2;

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

    UP = [1, 0; 0, 0];
    DOWN = [0, 0; 0, 1];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 针对六角晶格 Honeycomb，并考虑开边界为 y 方向，x 方向为周期结构
    % 开边界端原胞为 A-A 锯齿结构
    % 给定沿 y 方向有 n 个 Honeycomb 结构，该方向将包含 A-B 结构数 N = 2n

    N = n;

    H = zeros(4*N*2, 4*N*2);

    % 第一项（Hopping）

    H0 = zeros(4*N*2, 4*N*2);
    H0_mid = kron(eye(N),kron(sigmaz,kron(sigma0,[0, t*(exp(1i*(dot(k,a2))) + exp(1i*(dot(k,a3)))); 0, 0])));
    H0_hopping = circshift(eye(N),-1);
    H0_hopping(N,1) = 0;
    H0_hopping = kron(H0_hopping,kron(sigmaz,kron(sigma0,[0, 0; t*exp(-1i*(dot(k,a1))), 0])));
    H0 = H0_mid + H0_hopping;
    H0 = H0 + H0';



    % 第二项（SO）

    H_SO = zeros(4*N*2, 4*N*2);

    H_SO_a = kron(sigma0,kron(sigmaz,sigmaz.*(exp(1i*dot(k,b1))-exp(-1i*dot(k,b1)))));
    H_SO_mid = kron(eye(N),H_SO_a);
    H_SO_up = circshift(eye(N),-1);
    H_SO_up(N,1) = 0;
    H_SO_up = kron(H_SO_up,kron(sigma0,kron(sigmaz,sigmaz.*(exp(1i*dot(k,b2))-exp(-1i*dot(k,b3))))));
    H_SO_down = -H_SO_up';

    H_SO = -1i*lambda_SO*(H_SO_mid + H_SO_up + H_SO_down);


    % 第三项（mu）

    H_mu = zeros(4*N*2, 4*N*2);

    H_mu = kron(eye(N),kron(sigmaz,kron(sigma0,kron(sigma0, mu))));


    % 第四项(SC)
    H_SC0 = zeros(4*N*2,4*N*2);
    H_SC0_a = kron(sigmay,kron(sigmay,sigma0));
    H_SC0 = Delta0 * kron(eye(N),H_SC0_a);


    H_SC1 = zeros(4*N*2,4*N*2);

    H_SC1_mid = kron(eye(N),kron(sigmay,kron(sigmay,[0, t*(exp(1i*(dot(k,a2))) + exp(1i*(dot(k,a3)))); 0, 0])));
    H_SC1_hopping = circshift(eye(N),-1);
    H_SC1_hopping(N,1) = 0;
    H_SC1_hopping = kron(H_SC1_hopping,kron(sigmay,kron(sigmay,[0, 0; t*exp(-1i*(dot(k,a1))), 0])));
    H_SC1 = H_SC1_mid + H_SC1_hopping;
    H_SC1 = Delta1 * (H_SC1 + H_SC1');

    H_SC2 = zeros(4*N*2,4*N*2);

    H_SC2_up = circshift(eye(N),-1);
    H_SC2_up(N,1) = 0;
    H_SC2_a = kron(sigmay,kron(sigmay,sigma0));
    H_SC2_b1 = H_SC2_a(1:4,5:8).*( exp( 1i*(dot(k,b2))) + exp(-1i*dot(k,b3)) );
    H_SC2_b = [zeros(4,4), H_SC2_b1; zeros(4,4), zeros(4,4)];
    H_SC2_up = kron(H_SC2_up, H_SC2_b);
    H_SC2_down = H_SC2_up';
    H_SC2_mid = H_SC2_a(1:4,5:8).*exp(1i*dot(k,b1));
    H_SC2_mid = [zeros(4,4), H_SC2_mid; H_SC2_mid', zeros(4,4)];
    H_SC2_mid = kron(eye(N), H_SC2_mid);

    H_SC2 = Delta2 * (H_SC2_up + H_SC2_down + H_SC2_mid);


    H_SC = H_SC0 + H_SC1 + H_SC2;

    % 总的 Hamiltonian

    H = H0 + H_SO - H_mu - H_SC;

    % 以上是针对第一类边界，满足 A 开头，B 结尾
    % 接下来采用第二类边界，满足 A 开头，A 结尾
    % 等效于删除最后一个子块里与 B 相关的项，若不要请注释



    if type == 3
        H(4*(N-1)*2+8,:) = [];
        H(4*(N-1)*2+6,:) = [];
        H(4*(N-1)*2+4,:) = [];
        H(4*(N-1)*2+2,:) = [];
        H(:,4*(N-1)*2+8) = [];
        H(:,4*(N-1)*2+6) = [];
        H(:,4*(N-1)*2+4) = [];
        H(:,4*(N-1)*2+2) = [];
        H(7,:) = [];
        H(5,:) = [];
        H(3,:) = [];
        H(1,:) = [];
        H(:,7) = [];
        H(:,5) = [];
        H(:,3) = [];
        H(:,1) = [];
    elseif type == 4
        H(4*(N-1)*2+8,:) = [];
        H(4*(N-1)*2+6,:) = [];
        H(4*(N-1)*2+4,:) = [];
        H(4*(N-1)*2+2,:) = [];
        H(:,4*(N-1)*2+8) = [];
        H(:,4*(N-1)*2+6) = [];
        H(:,4*(N-1)*2+4) = [];
        H(:,4*(N-1)*2+2) = [];
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

    H = zeros(8*N, 8*N);

    % 第一项（Hopping）

    H0 = zeros(4*N, 4*N);
    H0_mid  = kron(eye(N), kron(sigmaz,kron(sigma0,[0, t*exp( 1i*(dot(k,a1))); 0, 0])));
    H0_up   = circshift(eye(N),-1);
    H0_up(N,1) = 0;
    H0_down = H0_up';
    H0_up   = kron(H0_up,  kron(sigmaz,kron(sigma0,[0, t*exp( 1i*(dot(k,a3))); 0, 0])));
    H0_down = kron(H0_down,kron(sigmaz,kron(sigma0,[0, t*exp( 1i*(dot(k,a2))); 0, 0])));
    H0 = H0_mid + H0_up + H0_down;
    H0 = H0 + H0';

    % 第二项（SO）

    H_SO = zeros(8*N, 8*N);
    H_SO_up = circshift(eye(N),-1);
    H_SO_up(N,1) = 0;
    H_SO_upup = circshift(H_SO_up,-1);
    H_SO_upup(N,2) = 0;

    H_SO_up = kron(H_SO_up,kron(sigmaz,kron(sigmaz,sigmaz*(exp(1i*dot(k,b2))+exp(1i*dot(k,b3))))));
    H_SO_down = -H_SO_up';
    H_SO_upup = kron(H_SO_upup,kron(sigmaz,kron(sigmaz,sigmaz*(-exp(-1i*dot(k,b1))))));
    H_SO_downdown = -H_SO_upup';

    H_SO = 1i*lambda_SO*(H_SO_up + H_SO_down + H_SO_upup + H_SO_downdown);

    % 第三项（mu）

    H_mu = zeros(4*N*2, 4*N*2);
    H_mu = kron(eye(N),kron(sigmaz,kron(sigma0,kron(sigma0, mu))));
    
    % 第四项（SC）

    H_SC0 = zeros(4*N*2,4*N*2);
    H_SC0 = zeros(4*N*2,4*N*2);

    H_SC0_a = -kron(sigmay,kron(sigmay,sigma0));
    H_SC0 = Delta0 * kron(eye(N),H_SC0_a);
    

    H_SC1_mid = -kron(eye(N),kron(sigmay,kron(sigmay,[0, exp(1i*(dot(k,a1))); 0, 0])));
    H_SC1_up = circshift(eye(N),-1);
    H_SC1_up(N,1) = 0;
    H_SC1_down = H_SC1_up';
    H_SC1_up   = -kron(H_SC1_up,  kron(sigmay,kron(sigmay,[0, exp(1i*(dot(k,a2))); 0, 0])));
    H_SC1_down = -kron(H_SC1_down,kron(sigmay,kron(sigmay,[0, exp(1i*(dot(k,a3))); 0, 0])));
    H_SC1 = H_SC1_mid + H_SC1_up + H_SC1_down;
    H_SC1 = Delta1 * (H_SC1 + H_SC1');

    H_SC2 = zeros(4*N*2,4*N*2);
    H_SC2_up = circshift(eye(N),-1);
    H_SC2_up(N,1) = 0;
    H_SC2_down = H_SC2_up';
    H_SC2_upup = circshift(H_SC2_up,-1);
    H_SC2_upup(N,2) = 0;
    H_SC2_downdown = H_SC2_upup';

    H_SC2_up   = -kron(H_SC2_up,  kron(sigmay,kron(sigmay,[0, (exp( 1i*(dot(k,b2)))+exp( 1i*(dot(k,b3)))); 0, 0])));
    H_SC2_down = -kron(H_SC2_down,kron(sigmay,kron(sigmay,[0, (exp(-1i*(dot(k,b2)))+exp(-1i*(dot(k,b3)))); 0, 0])));
    H_SC2_upup     = -kron(H_SC2_upup,    kron(sigmay,kron(sigmay,[0, exp(-1i*(dot(k,b1))); 0, 0])));
    H_SC2_downdown = -kron(H_SC2_downdown,kron(sigmay,kron(sigmay,[0, exp( 1i*(dot(k,b1))); 0, 0])));
    H_SC2 = H_SC2_up + H_SC2_upup + H_SC2_down + H_SC2_downdown;
    H_SC2 = Delta2 * (H_SC2 + H_SC2');
    
    H_SC = H_SC0 + H_SC1 + H_SC2;

    H = H0 + H_SO - H_mu + H_SC;

end

end




