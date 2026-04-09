function [x_sol, k, corr] = Gauss_Seidel_solver(A, b, kmax,tol, x_0)
    N = rows(A);
    n = N^(1/3);
    n = round(n);
    k = 0; x = x_0;
    do
        corr = 0;

        for i = 1:N
            xi = b(i);

            for j = max(i-n^2,1):min(i+n^2,N)

                if (j != i)
                    xi = xi - A(i, j) * x(j);
                endif

            endfor
            xi = xi / A(i, i);
            corr += (xi - x(i))^2;
            x(i) = xi;
        endfor
        corr = sqrt(corr); k = k + 1;
    until (corr < tol || k > kmax)
    x_sol = x;
endfunction;
