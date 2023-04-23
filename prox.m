function x = prox(problem, x, varargin)
    switch problem
        case 'l0'
            x = (varargin{1} < 0.5*x.^2) .* (x);
        case 'l1'
            x = sign(x) .* max(abs(x)-varargin{1},0);
        case 'i+'
            x = max(x, 0);
        otherwise
            warning('Unexpected proximal problem!');
    end
end