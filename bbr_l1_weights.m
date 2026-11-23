function a = bbr_l1_weights(beta, tau, m)
%BBR_L1_WEIGHTS  L1 weights (2.3a): a(i+1) = a_i^{(beta)}, i = 0..m.
%   For beta = 1 this gives a_0 = 1/tau and a_i = 0 (i >= 1).

i = (0:m)';
if beta == 1
    a = zeros(m+1,1);
    a(1) = 1/tau;
else
    a = ((i+1).^(1-beta) - i.^(1-beta)) / (tau^beta * gamma(2-beta));
end
end
