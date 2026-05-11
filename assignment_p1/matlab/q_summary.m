% Print a single block with every numeric answer used in the report.
parameters;

fprintf('\n=========== KEY DERIVED QUANTITIES ===========\n');
fprintf('I_beam      = %.4e m^4\n', I_beam);
fprintf('k_per_beam  = %.4f N/m\n', k_beam);
fprintf('k (4 beams) = %.4f N/m\n', k);
fprintf('m_top       = %.4e kg\n', m_top);
fprintf('m_centre    = %.4e kg\n', m_centre);
fprintf('m_combs     = %.4e kg\n', m_combs);
fprintf('m_total     = %.4e kg\n', m);
fprintf('omega_0     = %.4e rad/s\n', w0);
fprintf('f_0         = %.2f kHz\n', f0/1e3);
fprintf('\n--- Q3: gravity ---\n');
xg = m*g_grav/k;
fprintf('x_g         = %.3e m  = %.3f pm\n', xg, xg*1e12);

fprintf('\n--- Q4: damping & settling ---\n');
fprintf('A_m         = %.3e m^2\n', Am);
fprintf('b_air       = %.3e N s/m\n', b_air);
fprintf('b_clmp      = %.3e N s/m\n', b_clmp);
fprintf('b_total     = %.3e N s/m\n', b_tot);
fprintf('Q_air       = %.2f\n', Q_air);
fprintf('Q_clmp      = %.2f\n', Qclmp);
fprintf('Q_total     = %.2f\n', Q_tot);
fprintf('zeta        = %.5f\n', zeta);
ts = -log(0.05)/(zeta*w0);
fprintf('t_s (5%%)    = %.3e s = %.2f us\n', ts, ts*1e6);

fprintf('\n--- Q5: max acceleration ---\n');
A_DC  = g0*w0^2;
A_res = g0*w0^2/Q_tot;
fprintf('A_max,DC    = %.3e m/s^2 = %.0f g\n', A_DC, A_DC/g_grav);
fprintf('A_max,res   = %.3e m/s^2 = %.0f g\n', A_res, A_res/g_grav);

fprintf('\n--- Q6: phase shift at V=50 V ---\n');
V = 50;
xDC = 0.5*V^2*dCdx/k;
g_op = g0 + xDC;
neff = @(gg) (neff0-neffinf).*exp(-bcoef.*gg) + neffinf;
dphi = (2*pi*L/lambda) * (neff(g0) - neff(g_op));
fprintf('dC/dx       = %.3e F/m\n', dCdx);
fprintf('x_DC(50 V)  = %.3e m = %.3f nm\n', xDC, xDC*1e9);
fprintf('g(50 V)     = %.3f nm\n', g_op*1e9);
fprintf('Dphi(50 V)  = %.4f rad\n', dphi);

fprintf('\n--- Q7: pull-in (with finger-tip parallel plate) ---\n');
gt0 = Lc - d0;                        % rest tip gap
% solve g_t0 = u/2 + u^3/(2*h*Wc) + ... wait derive freshly
% Force balance at pull-in: g_t0 = 3u/2 + u^3/(2*h*Wc), u=g_t0-x_pi
fun = @(u) gt0 - 1.5*u - u.^3/(2*h*Wc);
u_pi = fzero(fun, gt0/2);
x_pi = gt0 - u_pi;
Vpi  = sqrt( k*u_pi^3 / (N*eps0*Wc*t) );
fprintf('g_t0        = %.3e m\n', gt0);
fprintf('u_pi (gap)  = %.3e m\n', u_pi);
fprintf('x_pi        = %.3e m\n', x_pi);
fprintf('V_pi        = %.2f V\n', Vpi);

fprintf('\n--- Q8: sensitivity at V_DC=10 V ---\n');
VDC = 10;
eta = VDC * dCdx;
Sphi = (2*pi*L/lambda) * (neff0-neffinf) * bcoef * exp(-bcoef*g0);
S_dc  = Sphi * eta / k;
S_res = Sphi * eta * Q_tot / k;
fprintf('eta         = %.3e N/V\n', eta);
fprintf('S_phi       = %.3e rad/m\n', Sphi);
fprintf('|H|_DC      = %.3e rad/V = %.3f mrad/V\n', S_dc, S_dc*1e3);
fprintf('|H|_res     = %.3e rad/V = %.3f mrad/V\n', S_res, S_res*1e3);

fprintf('\n--- Q8 plot end-points: low-f sensitivity at 5 kV ---\n');
S_5k = Sphi * (5e3 * dCdx) / k;
fprintf('|H|_DC@5kV  = %.3e rad/V = %.3f rad/V\n', S_5k, S_5k);

fprintf('==============================================\n\n');
