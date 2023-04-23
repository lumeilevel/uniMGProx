function gnorm = gradmap(problemF, problemG, x, L, varargin)
    gnorm = L * norm(x - prox(problemG, x-grad(problemF, x, varargin)/L, 1/L));
    switch problem
        case 'l0'
            gnorm = (varargin{2} < 0.5*varargin{1}^2) .* (varargin{1});
        case 'l1'
            gnorm = sign(varargin{1}) .* max(abs(varargin{1})-varargin{2},0);
        case 'indicator'
            gnorm = max(varargin{1}, 0);
        otherwise
            warning('Unexpected problem!');
    end
end