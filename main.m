clear;
N = [50 64 82 100 128 150 200 256];
N2 = [100, 256, 512, 1024, 2048, 4096, 8142, 10000];
pro = struct('Fn', 'QP', 'Gn', 'i+', 'lambda', 1e-2);
par = struct('tol', 1e-15, 'verbose', 1);
smooth = [1, 5, 10];    lsm = length(smooth);
levels = [1];           llv = length(levels);
choice = [2];           lch = length(choice);
lex = lsm * llv + 1;
infObj = 1e5;

for k = 1 : lch
    if strcmp(pro.Fn, 'QP')
        n = N(choice(k));
        h = 1 / (n + 1);
        pro.A = gallery('poisson', n) / h^2;
        pro.L0 = 8 / h^2;
        phi = max(sin((1 : n)*3*pi/n), 0);
        phi = vec(phi' * phi);
        pro.b = pro.A * phi;
        par.x_ini = rand(n^2, 1);
    else
        n = N2(k);
        pro.A = randn(n);
        pro.b = randn(n,1);
        pro.L0 = norm(A'*A);
        par.x_ini = randn(n,1);
        par.tol = 1e-9;
    end
    if strcmp(pro.Fn, 'LRbias')
        pro.c = randn(n,1);
    else
        pro.c = zeros(n,1);
    end
    par.max_iter = n^2 * 1e1;
    mu0 = svds(pro.A, 1, 'smallest');
    conv_fact = 1 - mu0 / pro.L0;
%     L = floor(2*log2(n)) - 1;
    
    x = cell(lex, lch);
    hist = cell(lex, lch);
    disp('APG');
    tic;
    [x{1}, hist{1}] = apg(pro, par);
    toc;
    infObj = min(infObj, min(hist{1+(k-1)*lex}.F));
    
    for i = 1 : lsm
        s = smooth(i);
        for j = 1 : llv
            L = levels(j);
            fprintf('\nMGProx-%d with level %d\n', s, L);
            t0 = tic;
            [x{i*j+1}, hist{i*j+1}] = mgproxL(A, p0, L0, x_ini, tol, L, s, verbose);
            toc(t0);
            infObj = min(infObj, min(hist{i*j+1+(k-1)*lex}.F));
        end
    end
end

% tic;
% cvx_begin quiet
%     variable y(N^2)
%     minimize (0.5*quad_form(y,Q0) + p0'*y)
%     subject to
%         y >= 0;
% cvx_end
% toc;