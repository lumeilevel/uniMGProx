function gnorm = gradmap(problemF, problemG, x, L, varargin)
    gnorm = L * norm(x - prox(problemG, x-grad(problemF, x, varargin{1:end-1})/L, varargin{end}/L));
end