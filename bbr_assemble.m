function [K, M, area] = bbr_assemble(p, t)
%BBR_ASSEMBLE  P1 stiffness matrix K and consistent mass matrix M.

x = reshape(p(t,1), [], 3);
y = reshape(p(t,2), [], 3);
b = [y(:,2)-y(:,3), y(:,3)-y(:,1), y(:,1)-y(:,2)];
c = [x(:,3)-x(:,2), x(:,1)-x(:,3), x(:,2)-x(:,1)];
area = 0.5*(b(:,1).*c(:,2) - b(:,2).*c(:,1));

nT = size(t,1);
I = zeros(9*nT,1); J = I; Kv = I; Mv = I;
q = 0;
for r = 1:3
    for s = 1:3
        idx = q+1 : q+nT;
        I(idx)  = t(:,r);
        J(idx)  = t(:,s);
        Kv(idx) = (b(:,r).*b(:,s) + c(:,r).*c(:,s)) ./ (4*area);
        Mv(idx) = area/12 * (1 + (r==s));
        q = q + nT;
    end
end
N = size(p,1);
K = sparse(I, J, Kv, N, N);     % sparse() sums duplicate entries
M = sparse(I, J, Mv, N, N);
end
