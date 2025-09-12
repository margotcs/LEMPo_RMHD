% =========================================================================
% Example: Run LEMPo_ideal for m/n=1/1 mode with a reversed q profile
% =========================================================================
%
% This script demonstrates how to run the LEMPo_ideal solver using a 
% strongly reversed q profile as input. 
% This example includes shaping.
% We build the profile and opts struct for inputs to LEMPo_ideal, and as 
% outputs we get plots of the radial plasma displacement, and the eigenvalue 
% \gamma / \omega_A.
%
% -------------------------------------------------------------------------

clear all; 

%% ------------------- Mode and Grid Selection ----------------------------

m = 1;           % Poloidal mode number
n = 1;           % Toroidal mode number
N = 3000;        % Number of grid points (higher N for better resolution)
toPlot = 1;      % Whether to show output plots

%% ------------------- Options for the Solver -----------------------------
% you can leave the opts struct empty : opts = [], the code will still run.

opts.additionalPlot = 0;     % Show profile/grid plots before solving 
opts.addShaping = 1;         % Shaping of flux surfaces
opts.R0 = 3;                 % Major radius (default is R0=3)

%% ------------------- Strongly reversed q Profile Setup ------------------

x = sym('x') % x stands for the radial coordinate
% Defining the reversed q profile
qs = m/n;
q0 = 1.2;             % q at axis
rA = 0.655;           
muZero = 3.8824;      
mu1 = 0;              
f1 = -0.238;          
r11 = 0.4286;         
r12 = 0.304;          

mu = @(x) muZero + mu1.*x.^2;
F1 = @(x) 1 + f1.*exp(- ((x - r11)/r12).^2);
r0 = rA*abs((m/(n*q0))^(mu(rA)) - 1)^(-1/(2*mu(rA)));

profiles.q = @(x) q0.*F1(x).*(1 + (x/r0).^(2*mu(x))).^(1./mu(x));

%% ------------------- Pressure Profile Setup -----------------------------

RS = 0.2472;          % We force the value of alpha at the 1st rat. surf.
alphaWanted = 0.1;    % Desired alpha at RS
beta0 = 0.01;         % Central beta
R0 = opts.R0;         % Major radius

% Custom pressure profile (for a parabolic one, just use profiles.alpha = ..
beta1 = alphaWanted/(5*beta0*qs^2*R0);
profiles.beta = @(x) -beta0.*tanh(10.*beta1.*(x-RS))./2 + beta0/2;

%% ------------------- Shaping Parameters ---------------------------------
        
profiles.elong = 1.2;         % Elongation \kappa, considered constat
profiles.triangByR = -0.2;    % Triangularity/ r = \delta / r, considered constant


%% ------------------- Initial Eigenvalue Guess ---------------------------

% Initial guess for the eigenvalue, you should
% provide something much above the expected result.
ev_guess = 1e-3;      

%% ------------------- Run the Solver -------------------------------------

[ev, result] = LEMPo_ideal(m, n, N, ev_guess, toPlot, profiles, opts);

% -------------------------------------------------------------------------
% Outputs:
%   - ev: Growth rate (gamma / omega_A)
%   - result: Struct containing eigenvector(s), grid, profiles, etc.
% -------------------------------------------------------------------------

%% ------------------- End of Example -------------------------------------