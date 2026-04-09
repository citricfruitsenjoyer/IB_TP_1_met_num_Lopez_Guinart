fid = fopen("times_rnd_pivot.csv", "a");  % open file for writing
fid
rng(123);   % 123 es la semilla que quieras
fprintf(fid,'N, t, err,norm_x\n');
% A=rand(10,10);
fclose(fid);

for k=1:12;
    fid = fopen("times_rnd_pivot.csv", "a");  % open file for writing
    N = 2^k; 
    A=rand(N,N);
    x = rand(N,1);
    b = A*x;
    factor=10^-12;
    mask = rand(N,N) < 0.3;
    A(mask) = A(mask) * factor; % Mucho mas eficiente que loop abajo
    tic;
    x_sol = LU_solver(A,b);

    t = toc;
    err = norm(x-x_sol);
    x_norm = norm(x);
    fprintf("N = %d, time = %f seconds,err = %f,x_norm = %f \n", N, t, err,x_norm);        % print to console
    fprintf(fid, "%d,%.16f ,%.16f,%.16f\n", N, t, err,x_norm);
    fclose(fid);
    % fprintf(fid, "N = %d, time = %.16f seconds, Err = \n", N, t,err);  % write to file
end

