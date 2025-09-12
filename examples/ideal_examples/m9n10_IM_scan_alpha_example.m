% =========================================================================
% Example: Run IM_LEMPo_ideal for a m/n=9/10 mode (scan on pressure)
% =========================================================================
%
% This script demonstrates how to run the IM_LEMPo_ideal solver for a 
% range of pressure gradients (alpha values) for a m/n=9/10 interchange mode.
% It plots the growth rates from IM_LEMPo_ideal and the theoretical Mercier 
% marginal stability limit.
%
% -------------------------------------------------------------------------

clear all;

%% ------------------- Mode and Grid Selection ----------------------------

m = 9;           % Poloidal mode number
n = 10;          % Toroidal mode number
N = 3000;        % Number of grid points (higher N for better resolution)
toPlot = 1;      % Whether to show output plots

%% ------------------- Options for the Solver -----------------------------

opts.additionalPlot = 0;    % Show profile/grid plots before solving
opts.N_tries = 15;          % Max number of tries for eigs solver

%% ------------------- Set up Scan and Profiles ---------------------------

N_it = 5;                                  % Number of iterations in the scan
alpha_vec = linspace(0.01, 0.2, N_it);     % Range of alpha values (pressure gradients)

% Using the default q profile (see in utils/getProfiles.m). For a custom q
% profile use a function handle profiles.q = @(x) ...


%% ---------------- Obtain Mercier Marginal Stability ---------------------


RS = 0.2457;           % Main rational surface location
s = 0.055246;          % Shear at rational surface (recompute if q profile changes)
iar_rs = RS / 3;       % Inverse aspect ratio at rational surface (R0 = 3)
qs = m / n;            % Safety factor at rational surface

% Calculate critical alpha for Mercier criterion
alpha_crit = iar_rs^(-1) * (s^2 / 4) * (qs^(-2) - 1)^(-1);

%% ------------------- Initial Eigenvalue Guess ---------------------------

% Initial guess for the eigenvalue, you should
% provide something much above the expected result.
ev_guess = 5e-3;       

%% ------------------- Run the Alpha Scan ---------------------------------

ev = zeros(N_it, 1);    

for i = 1:N_it

    profiles.alpha = alpha_vec(i);   % This will set a parabolic pressure 
    % profile, with desired alpha value at the rational. For a different
    % type of pressure profile, use a function handle profiles.beta = @(x) ..

    [ev(i), result] = IM_LEMPo_ideal(m, n, N, ev_guess, toPlot, profiles, opts);

end

%% ------------------- Plot Results ---------------------------------------

set(0, 'defaultTextInterpreter', 'tex');
set(0, 'DefaultLineLineWidth', 3);
set(groot, 'defaultLineMarkerSize', 15);
set(gca, 'TickLabelInterpreter', 'tex');
set(0, 'DefaultAxesFontSize', 25)

cl = EPFLcolors();   % EPFL color palette (in utils)

figure
plot(alpha_vec, ev, 'o--', 'Color', cl.canard);             
hold on
xline(alpha_crit, 'k--', 'LineWidth', 2)                     

legend('LEMPo_{ideal}^{IM}', 'Mercier marginal stability')
xlabel('\alpha')
ylabel('\gamma / \omega_A')

% -------------------------------------------------------------------------
% Outputs:
%   - ev: Growth rates (gamma / omega_A) for each alpha value
%   - Plots of the radial plasma displacement
% -------------------------------------------------------------------------

%% ------------------- End of Example -------------------------------------