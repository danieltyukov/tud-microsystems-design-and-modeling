% Q8 last bullet: low-frequency sensitivity vs V_DC, 0..5 kV.
parameters;
outdir = fullfile(fileparts(mfilename('fullpath')),'..','report','figures');

VDC = linspace(0, 5e3, 1001);
% Operating gap shifts with V_DC (still using small-signal for AC part):
xDC = 0.5*VDC.^2*dCdx/k;
g_op = g0 + xDC;

% sensitivity uses derivative evaluated at operating point
Sphi_op = (2*pi*L/lambda) * (neff0 - neffinf) * bcoef .* exp(-bcoef.*g_op);
eta_op  = VDC * dCdx;
H_lf    = Sphi_op .* eta_op / k;

% Pure linear part (ignoring operating-point shift): proportional to V_DC
H_lin   = ((2*pi*L/lambda)*(neff0-neffinf)*bcoef*exp(-bcoef*g0)) * (VDC*dCdx) / k;

% Practical limits for context
% -- mass touches waveguide when x_DC = g0 (worst-case static)
V_contact = sqrt( 2*k*g0/dCdx );
% -- pull-in voltage from q_summary (re-derive)
gt0 = Lc - d0;
fun = @(u) gt0 - 1.5*u - u.^3/(2*h*Wc);
u_pi = fzero(fun, gt0/2);
Vpi = sqrt( k*u_pi^3 / (N*eps0*Wc*t) );

figure('Color','w','Position',[100 100 760 460]);
plot(VDC, H_lin, '--','LineWidth',1.4,'Color',[0.5 0.5 0.5]); hold on;
plot(VDC, H_lf,  '-', 'LineWidth',1.8,'Color',[0 0.45 0.74]);
grid on; box on;
xline(V_contact,':','Mass-vs-waveguide contact', ...
    'Color',[0.85 0.33 0.10],'LabelHorizontalAlignment','center','FontSize',9);
xline(Vpi,':',sprintf('Comb-tip pull-in (%.0f V)',Vpi), ...
    'Color',[0.49 0.18 0.56],'LabelHorizontalAlignment','center','FontSize',9);
xlabel('V_{DC} (V)');
ylabel('|\Delta\phi / V_{AC}| at f \rightarrow 0  (rad/V)');
legend('Linear (V_{DC}\cdot dC/dx) / k', ...
       'With operating-point shift', ...
       'Location','northwest');
title('Low-frequency modulation sensitivity vs DC bias');
xlim([0 5e3]);

set(gca,'FontSize',11);
exportgraphics(gca, fullfile(outdir,'q8_sens_vs_VDC.pdf'), 'ContentType','vector');
fprintf('saved q8_sens_vs_VDC.pdf  (V_contact=%.0f V, V_pi=%.0f V)\n', V_contact, Vpi);
