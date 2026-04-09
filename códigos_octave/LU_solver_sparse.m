function x_sol = LU_solver_sparse(A, b)
    N = rows(A);
    n = N^(1/3);
    n = round(n);

    U = A; % copiamos A en U, guardaremos en U las matrices intermedias A(k)
    L = eye(N,N); % L comienza siendo la matriz identidad deN ×N

    for k = 1:N - 1

        for i = k + 1: min(N, k+n^2)
            L(i, k) = U(i, k) / U(k, k); % calcula el coeficiente de L−1
            % opera sobre la fila i, recorriendo con j la parte queNo tiene ceros
            for j = k + 1:min(N, k+n^2)
                U(i, j) = U(i, j) - L(i, k) * U(k, j);
            endfor

            U(i, k) = 0; % simple "prolijidad"
        endfor

    endfor

    for k = 1:N
        y(k) = b(k);
        for i = k + 1:N
            b(i) = b(i) - L(i, k) * y(k); % pasa al 2do miembro el t´ermino de y(k)
        endfor

    endfor

    for k = N:-1:1
        x(k) = y(k) / U(k, k);
        for i = 1:k - 1
            y(i) = y(i) - U(i, k) * x(k); % pasa al 2do miembro el t´ermino de y(k)
        endfor
    endfor

    % toc
    x_sol = x;


endfunction;