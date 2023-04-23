function F = obj(problemF, problemG, x, varargin)
    % ends at varargin{2} or varargin{3}
    switch problemF
        case 'LR'
            F = 0.5*(varargin{1}*x-varargin{2})'*(varargin{1}*x-varargin{2});
        case 'QP'
            F = 0.5 * x' * varargin{1} * x + varargin{2}' * x;
        case 'LRbias'
            F = 0.5*(varargin{1}*x-varargin{2})'*(varargin{1}*x-varargin{2}) + varargin{3}'*x;
        otherwise
            warning('Unexpected smooth term!');
    end
    % starts from varargin{4}
    switch problemG
        case 'l0'
            F = F + varargin{4} * sum(x > eps);
        case 'l1'
            F = F + varargin{4} * norm(x, 1);
        case 'i+'
            
        otherwise
            warning('Unexpected nonsmooth term!');
    end
end