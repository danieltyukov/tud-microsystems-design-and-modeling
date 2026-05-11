% Q4: visualize the step response and the 5% settling band.
parameters;
outdir = fullfile(fileparts(mfilename('fullpath')),'..','report','figures');

% System: m xddot + b xdot + k x = F  =>  X(s)/F(s) = 1/(ms^2 + bs + k)
sys = tf(1, [m b_tot k]);
F0  = 1e-9;                   % arbitrary 1 nN step (linear scaling)
Tend = 5e-4;
t = linspace(0, Tend, 4000);
x = lsim(sys, F0*ones(size(t)), t);

x_inf = F0/k;
ts = -log(0.05)/(zeta*w0);

figure('Color','w','Position',[100 100 760 420]);
plot(t*1e6, x*1e9, 'LineWidth',1.6, 'Color', [0 0.45 0.74]); hold on;
yline(x_inf*1e9*1.05,'--','Color',[0.6 0.6 0.6]);
yline(x_inf*1e9*0.95,'--','Color',[0.6 0.6 0.6]);
yline(x_inf*1e9,':','Color',[0.4 0.4 0.4]);
xline(ts*1e6,':','Color',[0.85 0.33 0.10],'LineWidth',1.4, ...
    'Label',sprintf('t_s = %.0f \\mus',ts*1e6),'LabelVerticalAlignment','top');
xlabel('time (\mus)'); ylabel('displacement x_m (nm)');
title('Step response of suspended mass to a 1 nN step force');
grid on; box on;
xlim([0 Tend*1e6]);

set(gca,'FontSize',11);
exportgraphics(gca, fullfile(outdir,'q4_step_response.pdf'), 'ContentType','vector');
fprintf('saved q4_step_response.pdf\n');
