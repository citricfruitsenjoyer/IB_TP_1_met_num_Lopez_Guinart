function x_sol = LU_solver(A, b)
    % tic
    n = rows(A);
    p = 1:n;

    % for i = 1:n - 1
    %     [val, inx] = max(abs(A(i, i:end)));
    %     inx = inx + i - 1;
    %     A(:, [i inx]) = A(:, [inx i]);
    %     p([i inx]) = p([inx i]);
    % endfor
    % disp(p)
    U = A; % copiamos A en U, guardaremos en U las matrices intermedias A(k)
    L = eye(n, n); % L comienza siendo la matriz identidad de n × n

    for k = 1:n - 1
        [val, inx] = max(abs(U(k, k:end)));
        inx = inx + k - 1;
        U(:, [k inx]) = U(:, [inx k]);
        p([k inx]) = p([inx k]);
        for i = k + 1:n
            L(i, k) = U(i, k) / U(k, k); % calcula el coeficiente de L−1
            % opera sobre la fila i, recorriendo con j la parte que no tiene ceros
            for j = k + 1:n
                U(i, j) = U(i, j) - L(i, k) * U(k, j);
            endfor

            U(i, k) = 0; % simple "prolijidad"
        endfor

    endfor

    for k = 1:n
        y(k) = b(k);
        for i = k + 1:n
            b(i) = b(i) - L(i, k) * y(k); % pasa al 2do miembro el t´ermino de y(k)
        endfor

    endfor

    for k = n:-1:1
        x(k) = y(k) / U(k, k);
        for i = 1:k - 1
            y(i) = y(i) - U(i, k) * x(k); % pasa al 2do miembro el t´ermino de y(k)
        endfor

    endfor

    x_perm = x;
    x = zeros(n,1);
    x(p) = x_perm;
    % toc
    x_sol = x;
endfunction;