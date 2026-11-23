function out = bbr_solve(beta, lam, f, T, m, mesh, opts)
%BBR_SOLVE  BBR-type scheme (2.12) with the L1 discretization in time.
%
%   out = bbr_solve(beta, lam, f, T, m, mesh, opts)
%
%   beta : order in (0,1] (beta = 1 gives the classical BBR scheme)
%   lam  : Yosida parameter lambda > 0
%   f    : function handle f(t, x, y), vectorized in x, y
%   T, m : final time and number of time steps (tau = T/m)
%   mesh : struct with fields p, t, bnd, K, M (from bbr_disk_mesh/bbr_assemble)
%   opts : optional struct
%          .Phi        nonlinearity (default: positive part)
%          .checkStep  check (2.15) (default: true)
%          .saveG      vector of step indices k at which g_k is stored (default: [])
%
%   One step k = 1..m:
%       (M + lam*K) v_k = M (u_{k-1} + lam*f_k)       (Dirichlet 0)
%       g_k = Phi((v_k - u_{k-1}) / lam)              (= delta_tau^beta u_k)
%       u_k = ( g_k + sum_{j=1}^{k-1} (a_{j-1}-a_j) u_{k-j} ) / a_0      (2.3b)
%
%   out.U   : nI x (m+1), column k+1 = u_k on the interior nodes
%   out.G   : nI x numel(saveG), g_k at the requested steps
%   out.inner, out.tau, out.a

if nargin < 7, opts = struct(); end
if ~isfield(opts, 'Phi'),       opts.Phi = @(s) max(s, 0); end
if ~isfield(opts, 'checkStep'), opts.checkStep = true;     end
if ~isfield(opts, 'saveG'),     opts.saveG = [];           end

p = mesh.p; K = mesh.K; M = mesh.M;
N = size(p,1);
inner = setdiff((1:N)', mesh.bnd);
nI = numel(inner);
tau = T/m;

if opts.checkStep
    tm = bbr_step_condition(beta, lam);
    if tau > tm*(1 + 1e-12)
        error('tau = %.3e violates (2.15): need tau <= %.3e', tau, tm);
    end
end

a = bbr_l1_weights(beta, tau, m);      % a(1) = a_0, ..., a(m+1) = a_m
bcoef = a(1:end-1) - a(2:end);         % bcoef(j) = a_{j-1} - a_j, j = 1..m

% factorize M + lam*K once (sparse Cholesky with AMD ordering)
A  = M(inner,inner) + lam*K(inner,inner);
pa = amd(A);
R  = chol(A(pa,pa));
MI = M(inner,:);

% 4-point Gauss-Legendre rule on [-1,1] for f_k = (1/tau) int f
gx = [-0.861136311594053; -0.339981043584856; 0.339981043584856; 0.861136311594053];
gw = [ 0.347854845137454;  0.652145154862546; 0.652145154862546; 0.347854845137454];

U = zeros(nI, m+1);                    % U(:,1) = u_0 = 0
G = zeros(nI, numel(opts.saveG));
ufull = zeros(N,1);
for k = 1:m
    t0 = (k-1)*tau; tc = t0 + tau/2;
    fk = zeros(N,1);
    for q = 1:4
        fk = fk + gw(q) * f(tc + 0.5*tau*gx(q), p(:,1), p(:,2));
    end
    fk = fk/2;

    ufull(:) = 0;
    ufull(inner) = U(:,k);
    rhs = MI*(ufull + lam*fk);
    v = zeros(nI,1);
    v(pa) = R \ (R' \ rhs(pa));

    g = opts.Phi((v - U(:,k))/lam);

    if k > 1
        hist = U(:,2:k) * bcoef(k-1:-1:1);   % sum_{i=1}^{k-1} (a_{k-i-1}-a_{k-i}) u_i
    else
        hist = 0;
    end
    U(:,k+1) = (g + hist) / a(1);

    js = find(opts.saveG == k);
    if ~isempty(js), G(:,js) = repmat(g, 1, numel(js)); end
end

out = struct('U', U, 'G', G, 'inner', inner, 'tau', tau, 'a', a);
end
