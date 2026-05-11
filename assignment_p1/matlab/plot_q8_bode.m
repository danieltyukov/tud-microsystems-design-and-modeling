% Q8: Bode diagram of Dphi(s)/V_AC(s) at V_DC = 10 V.
parameters;
outdir = fullfile(fileparts(mfilename('fullpath')),'..','report','figures');

VDC  = 10;
eta  = VDC * dCdx;                                 % electromech coupling [N/V]
g_op = g0 + 0.5*VDC^2*dCdx/k;                      % operating gap (~= g0)
Sphi = (2*pi*L/lambda) * (neff0-neffinf) * bcoef * exp(-bcoef*g_op);

% Transfer function: H(s) = Sphi * eta / (m s^2 + b s + k)
sys = tf( Sphi*eta, [m b_tot k] );

f = logspace(2, 8, 2000);          % 100 Hz .. 100 MHz
w = 2*pi*f;
[mag, phs] = bode(sys, w);
mag = squeeze(mag);
phs = squeeze(phs);

figure('Color','w','Position',[100 100 800 600]);

subplot(2,1,1);
loglog(f, mag, 'LineWidth', 1.8, 'Color', [0 0.45 0.74]);
hold on; grid on; box on;
S_dc  = Sphi*eta/k;
S_res = Sphi*eta*Q_tot/k;
yline(S_dc, '--', sprintf('|H|_{DC} = %.3g rad/V', S_dc), ...
    'LabelHorizontalAlignment','left','Color',[0.5 0.5 0.5]);
plot(f0, S_res, 'ro','MarkerFaceColor','r','MarkerSize',7);
text(f0, S_res, sprintf('  resonance: %.3g rad/V @ %.1f kHz', S_res, f0/1e3),...
    'VerticalAlignment','middle','FontSize',10);
xlabel('Frequency (Hz)');
ylabel('|\Delta\phi / V_{AC}|  (rad/V)');
title(sprintf('Modulation transfer function magnitude  (V_{DC} = %g V)',VDC));
xlim([f(1) f(end)]);

subplot(2,1,2);
semilogx(f, phs, 'LineWidth', 1.8, 'Color', [0.85 0.33 0.10]);
grid on; box on;
xlabel('Frequency (Hz)');
ylabel('\angle(\Delta\phi / V_{AC})  (deg)');
xlim([f(1) f(end)]);
yticks(-180:45:0);

exportgraphics(gcf, fullfile(outdir,'q8_bode.pdf'), 'ContentType','vector');
fprintf('saved q8_bode.pdf\n');
