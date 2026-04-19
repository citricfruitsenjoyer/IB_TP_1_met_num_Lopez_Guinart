fid = fopen("times_rnd_no_pivot_factor.csv", "a");  % open file for writing
fid
rng(123);   % 123 es la semilla que quieras
fprintf(fid,'N,k,t,err,norm_x\n');
% A=rand(10,10);
fclose(fid);

for k=0:12;
    fid = fopen("times_rnd_no_pivot_factor.csv", "a");  % open file for writing
    N = 200; 
    A=rand(N,N);
    x = rand(N,1);
    factor=10^(-k);
    mask = rand(N,N) < 0.3;
    A(mask) = A(mask) * factor; % Mucho mas eficiente que loop abajo
    b = A*x;
    tic;
    x_sol = LU_solver_unoptimized(A,b);

    t = toc;
    err = norm(x-x_sol);
    x_norm = norm(x);
    fprintf("N = %d, k = %d, time = %f seconds,err = %f,x_norm = %f \n", N,k, t, err,x_norm);        % print to console
    fprintf(fid, "%d,%d,%.16f ,%.16f,%.16f\n", N, k , t, err,x_norm);
    fclose(fid);
    % fprintf(fid, "N = %d, time = %.16f seconds, Err = \n", N, t,err);  % write to file
end

