% =========================================================================
% Example: Run LEMPo_res for a m/n=9/10 infernal mode (scan on Lundquist
% number)
% =========================================================================
%
% This script demonstrates how to run the LEMPo_res solver for a scan in
% the Lundquist number (SL) for a m/n=9/10 infernal mode.
% It plots the linear growth rates from LEMPo_res, and the scaling with
% Lundquist number. We are near ideal marginal stability with alpha = 0.05,
% so we obtain gamma / omega_A \propto SL^{-0.39} (approximately). This is
% slightly different from the resistive interchange scaling (SL^(-1/3))
% because we use a flat q profile and pick up infernal modes.
%
% -------------------------------------------------------------------------

clear all;

%% ------------------- Mode and Grid Selection ----------------------------

m = 9;              % Poloidal mode number
n = 10;             % Toroidal mode number
N = 3000;           % Number of grid points (higher N for better resolution)
toPlot = 1;         % Whether to show output plots

%% ------------------- Options for the Solver -----------------------------

opts.additionalPlot = 0;    % Show profile/grid plots before solving
opts.N_tries = 15;          % Max number of tries for eigs solver

%% ------------------- q Profile Setup ------------------------------------

qs = m/n;           % Safety factor at rational surface
q0 = 0.895;         % q at axis
rs = 0.3;           % End of low shear zone (location of RS)
d = 9.3;            % Parameter to obtain s ≈ 0.1 at rational surface
c = ((qs/q0)^d - 1) / rs^(2*d);
profiles.q = @(x) q0 * (1 + c * x.^(2*d)).^(1/d);

%% ------------------- Pressure/Alpha Profile -----------------------------

profiles.alpha = 0.05;       % This will set a parabolic pressure 
% profile, with desired alpha value at the rational. For a different
% type of pressure profile, use a function handle profiles.beta = @(x) ...

% alpha = 0.05 corresponds rougly to ideal marginal stability

%% ------------------- SL Scan Parameters ---------------------------------

N_it = 5;                          
SL = logspace(6, 12, N_it);         % Lundquist number values

%% ------------------- Initial Eigenvalue Guess ---------------------------

% Initial guess for the eigenvalue, you should
% provide something much above the expected result.
ev_guess = 5e-3;                     


%% ------------------- Run the SL Scan ------------------------------------

ev = zeros(N_it, 1);       

for i = 1:N_it
    profiles.SL = SL(i);             
    [ev(i), result] = LEMPo_res(m, n, N, ev_guess, toPlot, profiles, opts);
end

%% --------------- Save Results, and Fit ----------------------------------

save('ev_LEMPo_res.mat', 'SL', 'ev');


Const = polyfit(log(SL),log(ev), 1);
m = Const(1)
k = Const(2);
evFit = SL.^m.*exp(k);



%% ------------------- Plot Results ---------------------------------------

set(0, 'defaultTextInterpreter', 'tex');
set(0, 'DefaultLineLineWidth', 3);
set(groot, 'defaultLineMarkerSize', 15);
set(gca, 'TickLabelInterpreter', 'tex');
set(0, 'DefaultAxesFontSize', 25)

cl = EPFLcolors();   % EPFL color palette (in utils)

figure
loglog(SL, ev, 'o--', 'Color', cl.canard);                                
hold on
loglog(SL, evFit, ':', 'Color', cl.zinzolin);
legend('LEMPo_{res}', ['\proto S_L^{',num2str(m), '}'])
xlabel('S_L')
ylabel('\gamma / \omega_A')
title('Scaling of growth rate with Lundquist number S_L, for a m/n=9/10 mode')

% -------------------------------------------------------------------------
% Outputs of LEMPo_res:
%   - ev: Growth rates (gamma / omega_A) for each S_L value
%   - Plots of the radial plasma displacement
% -------------------------------------------------------------------------

%% ------------------- End of Example -------------------------------------