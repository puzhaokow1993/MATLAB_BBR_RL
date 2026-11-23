function taumax = bbr_step_condition(beta, lam, L)
%BBR_STEP_CONDITION  Largest tau allowed by (2.15):
%   Gamma(2-beta) tau^beta / (2(2 - 2^(1-beta))) <= lam/L.

if nargin < 3, L = 1; end
if beta == 1
    taumax = 2*lam/L;
else
    taumax = (lam/L * 2*(2 - 2^(1-beta)) / gamma(2-beta))^(1/beta);
end
end
