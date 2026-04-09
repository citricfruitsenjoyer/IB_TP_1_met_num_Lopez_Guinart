fid = fopen("times_GS_cubo.csv", "a");  % open file for writing
fid
rng(123);   % 123 es la semilla que quieras
fprintf(fid,'N,t,norm_x,k,err\n');
% A=rand(10,10);
fclose(fid);

for n=1:20;
    fid = fopen("times_GS_cubo.csv", "a");  % open file for writing
    n;
    N=n^3;
    SD=sparse(1:N,1:N,6*ones(N,1),N,N);
    if1=setdiff(1:N-1, n:n:n*(n^2-1));
    SU1=sparse(if1,if1+1,-1*ones(1,(n-1)*n^2),N,N);
    [i, j]=find([ones(n*(n-1),n); zeros(n,n)]);
    if2=i+(j-1)*n^2;
    SU2=sparse(if2,if2+n,-1*ones(1,n*(n-1)*n),N,N);
    if3=1:n^2*(n-1);
    SU3=sparse(if3,if3+n^2,-1*ones(1,n^2*(n-1)),N,N);
    A=SD+SU1+SU1'+SU2+SU2'+SU3+SU3'; %
    b = ones(1,N)/(n + 1)^2 ;

    x_0 = rand(N,1);
    tic;
    [x_sol, k, corr] = Gauss_Seidel_solver(A, b, 1000,1e-8, x_0);
    k = round(k);
    t = toc;
    % err = norm(x-x_sol);
    x_norm = norm(x_sol);
    fprintf("N = %d, time = %f seconds, x_norm= %f \n", N, t,x_norm);        % print to console
    fprintf(fid, "%d,%.16f,%.16f,%d,%.16f\n", N, t, x_norm, k, corr); 
    fclose(fid);
    % fprintf(fid, "N = %d, time = %.16f seconds, Err = \n", N, t,err);  % write to file
end

