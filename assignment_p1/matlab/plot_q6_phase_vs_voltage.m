% Q6: phase shift Dphi as a function of DC bias V_d in [-50, 50] V.
parameters;
outdir = fullfile(fileparts(mfilename('fullpath')),'..','report','figures');

V = linspace(-50, 50, 401);

% Comb-drive force is independent of x_m (perpendicular field only):
%   F = (1/2) V^2 * (N eps0 t / h)         [N]
% Spring restoring force:  F = k x_m  =>  x_m = F/k
xDC = 0.5 * V.^2 * dCdx / k;        % [m]
g   = g0 + xDC;                     % [m]

neff_g  = (neff0 - neffinf).*exp(-bcoef.*g)  + neffinf;
neff_g0 = (neff0 - neffinf)*exp(-bcoef*g0)   + neffinf;
dphi    = (2*pi*L/lambda) .* (neff_g0 - neff_g);   % [rad]

figure('Color','w','Position',[100 100 720 420]);
plot(V, dphi*1e3, 'LineWidth', 1.8, 'Color', [0 0.45 0.74]);
grid on; box on;
xlabel('DC bias V_d (V)');
ylabel('\Delta\phi (mrad)');
title('Phase shift as a function of DC bias');
xlim([-50 50]);

% annotate end value
[mx, ix] = max(abs(dphi));
text(V(ix), dphi(ix)*1e3*0.6, sprintf('|\\Delta\\phi(\\pm50 V)| = %.1f mrad', mx*1e3), ...
    'HorizontalAlignment','center','FontSize',10,'BackgroundColor',[1 1 1 0.7]);

set(gca,'FontSize',11);
exportgraphics(gca, fullfile(outdir,'q6_phase_vs_V.pdf'), 'ContentType','vector');
fprintf('saved q6_phase_vs_V.pdf\n');
