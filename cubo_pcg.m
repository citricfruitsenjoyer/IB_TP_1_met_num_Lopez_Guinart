n = 50;
N = n^3;
SD = sparse(1:N, 1:N, 6 * ones(N, 1), N, N);
if1 = setdiff(1:N - 1, n:n:n * (n^2 - 1));
SU1 = sparse(if1, if1 + 1, -1 * ones(1, (n - 1) * n^2), N, N);
[i, j] = find([ones(n * (n - 1), n); zeros(n, n)]);
if2 = i + (j - 1) * n^2;
SU2 = sparse(if2, if2 + n, -1 * ones(1, n * (n - 1) * n), N, N);
if3 = 1:n^2 * (n - 1);
SU3 = sparse(if3, if3 + n^2, -1 * ones(1, n^2 * (n - 1)), N, N);
A = SD + SU1 + SU1' + SU2 + SU2' + SU3 + SU3'; %

delta = 1 / (n + 1);
f = @(x, y, z) delta .* delta .* exp(-x) .* sin(pi .* y) .* sin(pi .* z);
b = [];

for k = 1:n

    for j = 1:n

        for i = 1:n

            b = [b f(i * delta, j * delta, k * delta)];

        endfor

    endfor

endfor

kmax=1000;
tic
Tsol = pcg(A, b', 1e-8, kmax);
toc
max(Tsol)



dlmwrite("cubo.csv", x_sol, ",", "precision", 12);


