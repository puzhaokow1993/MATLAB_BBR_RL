%RUN_MAIN  BBR-type scheme (2.12) on B_1 with Phi = (.)_+ for several beta.
%   Produces Figure 1 (snapshots, all beta) and Figure 2 (L2-norm in time).
clear; close all;

T   = 1.0;
lam = 0.1;
f   = @(t, x, y) 10*sin(2*pi*t) + 10*x;
betas = [0.3 0.5 0.7 0.9 1.0];

[p, t, bnd] = bbr_disk_mesh(32);              % 3169 nodes
[K, M]      = bbr_assemble(p, t);
mesh = struct('p', p, 't', t, 'bnd', bnd, 'K', K, 'M', M);
N = size(p,1);

res = struct();
for ib = 1:numel(betas)
    beta = betas(ib);
    tm = bbr_step_condition(beta, lam);
    m  = max(1000, ceil(T/tm));
    m  = 100*ceil(m/100);                     % so that t = 0.01*j are time levels
    saveIdx = (0:100)*(m/100);
    tic;
    out = bbr_solve(beta, lam, f, T, m, mesh);
    el = toc;

    Ufull = zeros(N, m+1);  Ufull(out.inner,:) = out.U;
    l2 = sqrt(sum(Ufull .* (M*Ufull), 1));

    res(ib).beta = beta; res(ib).m = m; res(ib).tm = tm;
    res(ib).l2 = l2; res(ib).Us = Ufull(:, saveIdx+1);
    fprintf('beta=%.1f  m=%d  tau_max=%.2e  time=%.1fs  max L2=%.3f  L2(T)=%.3f  min diff L2=%.2e\n', ...
            beta, m, tm, el, max(l2), l2(end), min(diff(l2)));
end
save('bbr_main.mat', 'res', 'p', 't', 'bnd', 'lam', 'T');

%% ---------------------------------------------------------------- figures
ramp  = hex2rgb({'#fcfcfb','#cde2fb','#9ec5f4','#6da7ec','#3987e5','#256abf','#184f95','#0d366b'});
cmap  = interp1(linspace(0,1,size(ramp,1)), ramp, linspace(0,1,256));
ordc  = hex2rgb({'#86b6ef','#5598e7','#2a78d6','#1c5cab','#0d366b'});   % beta = 0.3 ... 1

% Figure 1: snapshots for all beta (rows) at t = 0.25, 0.5, 0.75, 1 (columns)
times = [0.25 0.5 0.75 1.0];
nb    = numel(betas);
vmax  = max(cellfun(@(U) max(U(:)), {res.Us}));
fig1 = figure('Units','centimeters','Position',[2 2 16.5 20.5],'Color','w');
th = linspace(0, 2*pi, 400);
for i = 1:nb
    r = res(i);
    for j = 1:numel(times)
        ax = subplot(nb, numel(times), (i-1)*numel(times) + j);
        idx = round(times(j)*100) + 1;
        patch('Faces', t, 'Vertices', p, 'FaceVertexCData', r.Us(:, idx), ...
              'FaceColor', 'interp', 'EdgeColor', 'none');
        hold on
        plot(cos(th), sin(th), 'Color', [0.53 0.53 0.53], 'LineWidth', 0.5);
        axis equal off; colormap(ax, cmap); clim([0 vmax]);
        if i == 1, title(sprintf('$t$ = %g', times(j)), 'FontWeight', 'normal', 'Interpreter', 'latex'); end
        if j == 1
            text(-1.35, 0, sprintf('\\beta = %g', r.beta), 'Rotation', 90, ...
                 'HorizontalAlignment', 'center');
        end
    end
end
cb = colorbar('Position', [0.92 0.3 0.015 0.4]);   % uses the last axes' colormap and limits
ylabel(cb, '$u_{\tau}(t,x)$', 'Interpreter', 'latex');

% Figure 2: L2-norm in time
fig2 = figure('Units','centimeters','Position',[2 2 13 8],'Color','w'); hold on
for ib = 1:numel(betas)
    r = res(ib);  tt = linspace(0, T, r.m + 1);
    plot(tt, r.l2, 'Color', ordc(ib,:), 'LineWidth', 1.8, ...
         'DisplayName', sprintf('\\beta = %g', r.beta));
    text(T + 0.01, r.l2(end), sprintf('\\beta = %g', r.beta), 'FontSize', 8);
end
xlabel('$t$', 'Interpreter', 'latex'); ylabel('$||u_{\tau}(t,\cdot)||_{L^{2}(\Omega)}$', 'Interpreter', 'latex');
xlim([0 T]); ylim([0 inf]); box off; grid on; set(gca, 'XGrid', 'off');

%% ---------------------------------------------------------------- helpers
function rgb = hex2rgb(c)
rgb = zeros(numel(c), 3);
for k = 1:numel(c)
    h = c{k};
    rgb(k,:) = [hex2dec(h(2:3)), hex2dec(h(4:5)), hex2dec(h(6:7))]/255;
end
end

