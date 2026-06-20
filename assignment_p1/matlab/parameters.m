% =====================================================================
% MEMS optical phase shifter -- shared parameter file
% Loaded by every plotting/analysis script in this folder.
% =====================================================================

% ----- Geometry -----
L   = 100e-6;     % modulation section length [m]
g0  = 200e-9;     % modulation gap at rest    [m]
Lb  = 75e-6;      % beam length               [m]
Wb  = 5e-6;       % beam width                [m]
Lm  = 40e-6;      % central mass length       [m]
Wm  = 30e-6;      % mass section width        [m]
h   = 1e-6;       % comb-drive gap            [m]
d0  = 19e-6;      % comb overlap at rest      [m]
Lc  = 20e-6;      % comb finger length        [m]
Wc  = 1e-6;       % comb finger width         [m]
t   = 200e-9;     % suspended thickness       [m]
N   = 40;         % number of comb gaps       [-]

% ----- Material -----
Esi  = 170e9;     % Young's modulus  [Pa]
rho  = 2330;      % density          [kg/m^3]

% ----- Optics -----
lambda  = 1550e-9;  % wavelength            [m]
bcoef   = 0.3e6;    % refractive index decay [1/m]  (0.3/um)
neff0   = 2.6;      % n_eff at g=0
neffinf = 2.4;      % n_eff at g=infinity

% ----- Damping -----
cair  = 16.15;    % air damping coefficient [N s / m^3]
Qclmp = 500;      % clamping quality factor

% ----- Constants -----
eps0 = 8.854187817e-12;   % [F/m]
g_grav = 9.81;            % [m/s^2]

% ----- Beam stiffness (4 fixed-guided beams in parallel) -----
%   Each beam bends in the in-plane direction (across width Wb).
%   I = t * Wb^3 / 12   (Wb is the dimension along bending direction).
I_beam = t * Wb^3 / 12;
k_beam = 12 * Esi * I_beam / Lb^3;   % per beam, fixed-guided
N_beams = 4;
k = N_beams * k_beam;                % total stiffness

% ----- Effective mass -----
%   Figure-faithful moving structure (corrects the first submission, which
%   omitted the comb-spine plate and counted 40 fingers -> 2.33 pg):
%     - modulation bar      L  x Wm  (top, near waveguide)
%     - central column      Lm x Wm
%     - comb-spine plate    L  x Wm  (the lower "Mass" block; was MISSING)
%     - 20 moving fingers   Lc x Wc  (N=40 counts comb GAPS, not fingers)
N_fingers = 20;                          % number of MOVING comb fingers (= N/2)
m_bar    = rho * L  * Wm * t;            % top modulation bar
m_centre = rho * Lm * Wm * t;            % central column
m_spine  = rho * L  * Wm * t;            % comb-spine mass plate (added)
m_combs  = rho * N_fingers * Lc * Wc * t;% moving comb fingers
m = m_bar + m_centre + m_spine + m_combs;

% ----- Resonance -----
w0 = sqrt(k/m);
f0 = w0/(2*pi);

% ----- Comb-drive electromechanical constant per V_DC -----
dCdx = N * eps0 * t / h;             % [F/m] -- constant for ideal comb

% ----- Damping coefficients -----
Am = 2*L*Wm + Lm*Wm + N_fingers*Lc*Wc;  % top-side moving area (two plates + 20 fingers)
b_air  = cair * Am;
b_clmp = sqrt(k*m) / Qclmp;
b_tot  = b_air + b_clmp;
Q_air  = sqrt(k*m) / b_air;
Q_tot  = sqrt(k*m) / b_tot;
zeta   = b_tot / (2*sqrt(k*m));
