% =========================================================================
% Example: Run LEMPo_ideal for a m/n=8/10 infernal mode 
% =========================================================================
%
% This script demonstrates how run the LEMPo_ideal solver, using a
% piecewise, non resonant q profile as input. 
% We build the profile and opts struct for inputs of LEMPo_ideal, and as 
% outputs we get plots of the radial plasma displacement and returns the 
% eigenvalue \gamma / \omega_A.
%
% -------------------------------------------------------------------------

clear all; 

%% ------------------- Mode and Grid Selection ----------------------------

m = 8;            % Poloidal mode number (main harmonic)
n = 10;           % Toroidal mode number
N = 3000;         % Number of grid points (higher N for better resolution)
toPlot = 1;       % Whether to show output plots 

%% ------------------- Options for the Solver -----------------------------
% you can leave the opts struct empty : opts = [], the code will still run.

opts.additionalPlot = 1;  % Show profile/grid plots before solving 

%% ------------------- Pressure Profile Setup -----------------------------

profiles.alpha = 0.1;     % Value of alpha at r_s (rational surface)
% For other pressure profiles, you can use profiles.beta = @(x) ...

%% ------------------- q Profile Setup ------------------------------------
x = sym('x') % x stands for the radial coordinate

% q(x) is piecewise: flat for r<rb, quadratic for r>rb.
rb = 0.3;                 % End of low shear zone.
a = 1;                    % Plasma minor radius (normalized to 1)
qa = 10;                  % q(r=1)
delta_q = +0.001;         % Small offset for q(r<rb)

profiles.RS = rb;         % Need to give the code a value for r_s, because
% there is no actual rational surface in the plasma, and q is piecewise
% continuous. 

% Define q(x)
nonZeroPart = @(x) (m/n + delta_q) + (qa - (m/n + delta_q)).*(x - rb).^2./(a-rb)^2;
profiles.q = @(x) (m/n + delta_q).*(x <= rb) + nonZeroPart(x).*(x > rb);

% Compute derivatives for qp and qpp
nonZeroPartp = matlabFunction(diff(nonZeroPart(sym('x'))));
nonZeroPartpp = matlabFunction(diff(nonZeroPartp(sym('x'))));

profiles.qp = @(x) 0.*(x <= rb) + nonZeroPartp(x).*(x > rb);
profiles.qpp = @(x) 0.*(x <= rb) + nonZeroPartpp(x).*(x > rb);

%% ------------------- Sideband Rational Surface --------------------------

% Since q is piecewise continuous, need to provide the solver with value of
% upper sideband rational surface (where q = (m+1)/n )
RSp = vpasolve(nonZeroPart(sym('x')) == (m+1)/n , sym('x'), [0,1]);
profiles.RSp = double(RSp(2));  % 

%% ------------------- Initial Eigenvalue Guess ---------------------------

% Initial guess for the eigenvalue, you should
% provide something much above the expected result.
ev_guess = 1e-3;

%% ------------------- Run the Solver -------------------------------------

[ev, results] = LEMPo_ideal(m, n, N, ev_guess, toPlot, profiles, opts);

% -------------------------------------------------------------------------
% Outputs:
%   - ev: Growth rate (gamma / omega_A)
%   - results: Struct containing eigenvector(s), grid, profiles, etc.
% -------------------------------------------------------------------------

%% ------------------- End of Example -------------------------------------