function[H] = Hamiltonian(k)

A = 1;
B = 1;
M = 1;
Delta = 1;

sigma0 = [1,  0;   0,  1];
sigmax = [0,  1;   1,  0];
sigmay = [0, 1i; -1i,  0];
sigmaz = [1,  0;   0, -1];

kx = k(1);
ky = k(2);

H = (M-2*B*(2 - cos(kx) - cos(ky))).*kron(sigma0, sigmaz) ...
    + A*sin(kx).*kron(sigmax, sigmaz) ...
    + A*sin(ky).*kron(sigma0, sigmay) ...
    + Delta*(cos(kx) - cos(ky)).*kron(sigmax, sigmax);

end