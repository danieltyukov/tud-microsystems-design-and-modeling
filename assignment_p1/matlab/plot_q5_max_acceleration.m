% Q5: plot maximum tolerable acceleration vs frequency.
parameters;
outdir = fullfile(fileparts(mfilename('fullpath')),'..','report','figures');

f = logspace(0, 8, 2000);
w = 2*pi*f;

% |X/A| = m / |k - m w^2 + j w b|
H_xa = m ./ abs(k - m*w.^2 + 1j*w*b_tot);
A_max = g0 ./ H_xa;          % m/s^2

figure('Color','w','Position',[100 100 760 420]);
loglog(f, A_max/g_grav,'LineWidth',1.8,'Color',[0 0.45 0.74]); hold on;
grid on; box on;
% mark DC value and resonance
A_DC  = g0*w0^2;
A_res = g0*w0^2/Q_tot;
yline(A_DC/g_grav,'--','Color',[0.6 0.6 0.6], ...
    'Label',sprintf('A_{max,DC} = %.2g g',A_DC/g_grav));
plot(f0, A_res/g_grav, 'ro','MarkerFaceColor','r','MarkerSize',7);
text(f0*1.2, A_res/g_grav, sprintf('A_{max,res} = %.2g g  @  %.1f kHz', ...
    A_res/g_grav, f0/1e3), 'FontSize',10);
xlabel('Frequency (Hz)');
ylabel('Max tolerable acceleration (g)');
title('Maximum harmonic acceleration before mass touches waveguide');

set(gca,'FontSize',11);
exportgraphics(gca, fullfile(outdir,'q5_amax.pdf'), 'ContentType','vector');
fprintf('saved q5_amax.pdf  (A_DC=%.0f g, A_res=%.0f g)\n', A_DC/g_grav, A_res/g_grav);
