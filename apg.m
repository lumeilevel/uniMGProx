function [x, hist] = apg(pro, par)%apg(A, b, c, lambda, pro.L0, par.x_ini, par.tol, par.verbose)
    t = 1;  t0 = 1;
    eta = 0.5; tauk = pro.L0;
    x0 = par.x_ini; x = x0;
%     if strcmp(pro.Fn, 'LRbias')
%         objold = obj('LRbias', pro.Gn, x0, pro.A, pro.b, pro.c, pro.lambda);
%     else
%         objold = obj(pro.Fn, pro.Gn, x0, pro.A, pro.b, pro.lambda);
%     end
    hist.F = zeros(par.max_iter, 1);
    hist.G = zeros(par.max_iter, 1);
    hist.relDist = zeros(par.max_iter, 1);
    hist.relObjdiff = zeros(par.max_iter, 1);
    hist.G(1) = gradmap(pro.Fn, pro.Gn, x0, pro.L0, pro.A, pro.b, pro.c);
    for iter = 2 : par.max_iter
        hist.F(iter) = objv(pro.Fn, pro.Gn, x, pro.A, pro.b, pro.c, pro.lambda);
        hist.G(iter) = gradmap(pro.Fn, pro.Gn, x, pro.L0, pro.A, pro.b, pro.c);
        hist.relDist(iter) = norm(x-x0) / norm(x);
        hist.relObjdiff(iter) = abs(hist.F(iter) - hist.F(iter-1)) / max(hist.F(iter), 1);       
        % stopping criterion
        if hist.G(iter) / hist.G(1) <= par.tol
            hist.F = hist.F(1:iter);
            hist.G = hist.G(1:iter);
            hist.relDist = hist.relDist(1:iter);
            hist.relObjdiff = hist.relObjdiff(1:iter);
            if par.verbose
                fprintf('\n APG early stopping--iteration: %d\n', iter);
                fprintf('[c] proximal first-order optimality condition satisfied\n')
            end
            break
        end
        if iter > 10
            if max(hist.relDist(iter), 0.1*hist.relObjdiff(iter)) < par.tol
                if par.verbose
                    fprintf("\n APG Early Stopping--iteration: %d\n", iter);
                    fprintf('[a] relDist < %3.2e\n', par.tol);
                    fprintf("norm(X-Xold,'fro')/norm(X,'fro') = %f\n", hist.relDist(iter));
                end
                hist.F = hist.F(1:iter);
                hist.G = hist.G(1:iter);
                hist.relDist = hist.relDist(1:iter);
                hist.relObjdiff = hist.relObjdiff(1:iter);
                break
            end
            if max(0.5*hist.relDist(iter), 100*hist.relObjdiff(iter)) < par.tol
                if par.verbose
                    fprintf("\n APG Early Stopping--iteration: %d\n", iter);
                    fprintf('[b] relObjdiff < %3.2e\n', 0.01*par.tol);
                end
                hist.F = hist.F(1:iter);
                hist.G = hist.G(1:iter);
                hist.relDist = hist.relDist(1:iter);
                hist.relObjdiff = hist.relObjdiff(1:iter);
                break;
            end
        end
        y = x + (t0-1)/t*(x-x0);
        tau = eta * tauk;
        gradY = grad(pro.Fn, y, pro.A, pro.b, pro.c);
        fY = objv(pro.Fn, 'i+', y, pro.A, pro.b, pro.c);
        for j = 1 : 1e2
            s = prox(pro.Gn, y - gradY / tau, pro.lambda / tau);
            sy = s - y;
            if objv(pro.Fn, 'i+', s, pro.A, pro.b, pro.c) + eps <= fY + sy'*(0.5*tau*sy+gradY)
                tauk = tau;
                break;
            else
                tau = min(tau/eta, pro.L0);
            end
        end
        x0 = x;
        x = s;
        t0 = t; t = 0.5*(1+sqrt(1+4*t^2));
    end
end