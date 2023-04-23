function g = grad(problem, x, varargin)
    switch problem
        case 'LR'
            g = varargin{1}' * (varargin{1}*x-varargin{2});
        case 'QP'
            g = varargin{1} * x + varargin{2};
        case 'LRbias'
            g = varargin{1}' * (varargin{1}*x-varargin{2}) + varargin{3};
        otherwise
            warning('Unexpected smooth term!');
    end
end