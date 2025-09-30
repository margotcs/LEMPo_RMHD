function [ev, result]=IM_LEMPo_res(m, n, N, ev_guess, toPlot, profiles, opts)

% Run the interchange resistive model (so without infernal corrections)
%
% INPUTS:
%   m, n                - Poloidal & toroidal mode numbers of the main
%                         harmonic you want to look at
%   N                   - Number of points in the (radial) grid
%   ev_guess            - Initial guess for the eigenvalue (should be much above the
%                         expected value)
%   toPlot              - Boolean: 1 if you want outputs plots, 0 otherwise
%   profiles            - Struct containing all necessary inputs, see
%                         documentation
%   opts                - Struct containing options, can be set to []
%
% OUTPUTS:
%   ev      - Eigenvalue of the mode \gamma / \omega_A
%   result  - Struct containing displacements, sidebands, etc...

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%              Load Profiles, and Parameters you might modify
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

resistiveSimu = 1;

input = getProfiles(m, n, resistiveSimu, profiles, toPlot, opts);
[q, qp, qpp, beta, betap, betapp, alpha, alphap, R0, B0, Gamma, rStar, RS] = ...
    deal(input.q, input.qp, input.qpp, ...
    input.beta, input.betap, input.betapp, input.alpha, input.alphap, input.R0, ...
    input.B0, input.Gamma, input.rStar, input.RS);

[additionalPlot, addShaping, higherOscModes, reverse_qprofile, N_tries, SL] = ...
    deal(input.additionalPlot, input.addShaping, input.higherOscModes, ...
    input.reverse_qprofile,  input.N_tries, input.SL);

% Multiplication factors for the variables on the plots
factorChi = 10;
factorDeltaXi = 1;

if addShaping
    error(['You set opts.addShaping but shaping of flux surfaces is not ',...
        'included in the Interchange Models.'])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%      Set up grid (bunched around rational surface of main harmonic)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xall = gridBunching(RS, [], [], N, 1, opts); % (utils function)
dx = diff(xall);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%      Definitions of Some Quantities and First Outputs 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


s = @(x) x./q(x).*qp(x);      % shear
iar = @(x) x./R0;             % inverse aspect ratio
% below is radial derivative of Shafranov shift
deltap = @(x) q(x).^2./((x+1e-10).^3).*cumtrapz(x,x.^3./(q(x).^2.*R0) + x.^2.*alpha(x)./(q(x).^2));  

Dq = @(x) m./(n.*q(x)) - 1;    % this is the quantity qs ( 1/q - 1/qs )
Dqp = @(x) -(m.*qp(x))./(n.*q(x).*q(x));
Dqpp = @(x) -(m.*qpp(x))./(n.*q(x).^2) + (2*m.*qp(x).^2)./(n.*q(x).^3) ; 

% Other quantities that we need when including compression effects :
mu0 = pi*4e-7;      % ok for B in tesla
rho = @(x) 1./(1+(1.*x).^2).^2.*(1-x.^15) ;
mass_density = @(x) 8e-6.*rho(x) ;   % value in kg per m3 
alfven = @(x) B0./(sqrt(mu0.*mass_density(x)).*R0);
sound = @(x) sqrt( (Gamma.*beta(x))./(2.*mu0.*mass_density(x))  ).* ( B0./R0);

% First outputs:

result.rgrid = xall;
result.RS = RS;
result.alphaRS = alpha(RS);
result.s = s(RS);
result.q = q(xall);
result.deltap = deltap(xall);
result.SL = SL;
result.omegaS = sound(xall);
result.omegaA = alfven(xall);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%      Plot Profiles (if requested)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if additionalPlot
    makeFirstPlots(input, xall, N, 0); % (utils function)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%      Initialise & Define Coefficients in the Equations 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DLX_0D = []; DLX_1D = []; DLX_2D = [];
DRX_0D = []; DRX_1D = []; DRX_2D = [];
DRC_0D = []; DRC_1D = []; DRC_2D = [];
DRXG_0D = []; DRXG_1D = []; DRXG_2D = [];
OLC_0D = []; OLC_1D = []; OLC_2D = [];
ORX_0D = []; ORX_1D = []; ORX_2D = [];
ORC_0D = []; ORC_1D = []; ORC_2D = [];
DXG_LXG_0D = []; DXG_LXG_1D = []; DXG_LXG_2D = [];
DXG_RXG_0D = []; DXG_RXG_1D = []; DXG_RXG_2D = [];
DXG_RC_0D = []; DXG_RC_1D = []; DXG_RC_2D = [];

defineIMresOperators(); % nested function which defines operators 
%                         initialised above, model-specific

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%      Discretize Operators
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Simple 1st order finite difference discretisation

DLX      = discretize_S(DLX_0D, DLX_1D, DLX_2D);
DRX      = discretize_S(DRX_0D, DRX_1D, DRX_2D);
DRC      = discretize_S(DRC_0D, DRC_1D, DRC_2D);
DRXG     = discretize_S(DRXG_0D, DRXG_1D, DRXG_2D);

OLC      = discretize_S(OLC_0D, OLC_1D, OLC_2D);
ORC      = discretize_S(ORC_0D, ORC_1D, ORC_2D);
ORX      = discretize_S(ORX_0D, ORX_1D, ORX_2D);

DXG_LXG  = discretize_S(DXG_LXG_0D, DXG_LXG_1D, DXG_LXG_2D);
DXG_RXG  = discretize_S(DXG_RXG_0D, DXG_RXG_1D, DXG_RXG_2D);
DXG_RC   = discretize_S(DXG_RC_0D, DXG_RC_1D, DXG_RC_2D);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%      Boundary Conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
applyBCs(); % nested function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%      Build the Matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Here we build A and B to solve A v = gamma/omegaA B v

Z = sparse(N,N); S = speye(N);

  A = [Z,    S,    Z,     Z,    Z;
       DRX,  Z,   DRC,   DRXG,  Z;
       ORX,  Z,   ORC,    Z,    Z;
       Z,    Z,    Z,     Z,    S;
       Z,    Z,  DXG_RC, DXG_RXG, Z];
   
  B = [S,    Z,    Z,     Z,    Z;
       Z,  DLX,    Z,     Z,    Z;
       Z,    Z,   OLC,    Z,    Z;
       Z,    Z,    Z,     S,    Z;
       Z,    Z,    Z,     Z,  DXG_LXG];

% Check matrices A & B for NaN values
if any(isnan(A(:)))
    if any(isnan(B(:)))
        error('Matrices A and B contain NaN values.');
    else
        error('Matrix A contains NaN values');
    end
else
   if any(isnan(B(:)))
       error('Matrix B contains NaN values');
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%      Solve Eigenvalue Problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[V, ev] = tryEigs(A, B, ev_guess, N_tries, xall, N, rStar, reverse_qprofile, ...
    higherOscModes, additionalPlot); % (utils function)

if isempty(V) % if couldnt get any growth rates: stop the execution
    ev = 0;
    return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%      Recover Eigenvalue, Normalise Eigenvectors, Build Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if sum(V(:)) < 0
    V(:) = -V(:);
end

Xi  = V(1:N);    maxXi = max(Xi);
Xi  = Xi./maxXi;
Chi = V(2*N+1:3*N);  Chi = factorChi*Chi./maxXi;
XiG = V(3*N+1:4*N); XiG = factorDeltaXi*XiG./maxXi;

result.ev = ev;
result.Xi = Xi;
result.Chi = Chi./factorChi; % the multiplying factor is just to plot
result.XiG = XiG./factorDeltaXi;
result.guess = ev_guess;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%      Plot Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if toPlot
    makePlots(); % nested function
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%      Nested Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function defineIMresOperators()
    % How are the elements named in this subfunction :
    %   X_ZD : 
    %   - X corresponds to the side of the dispersion relation 
    % we are looking at. DLX marks the terms which will multiply gamma^2
    % (inertia), it will be in the B matrix, when solving B  Xi  ev =
    % A Xi. DRX will be in the A matrix (contains Mercier term, FLB etc)
    %   - Z is 1 if we are looking at the first derivative of the variable
    %   given by Y, 2 for second derivative, etc.
    
% -------------------------------------------------------------------------
    % Coefficients in the main mode governing equation (they have an additional 
    % factor (-r) wrt the paper or the documentation):
% -------------------------------------------------------------------------
  
    % Xi contribution to the main harmonic governing equation:
    DLX_0D = @(x) x - x./m.^2 + (2.*x.*m.^2)./n.^2 - 2.*x/n.^2  ; %%%
    DLX_1D = @(x) -3.*x.^2.*m.^(-2) - 2.*x.^2./n.^2  ;    % 
    DLX_2D = @(x) -x.^3.*m.^(-2) - 2.*x.^3./n.^2  ;    % vient 
           
    DRX_0D = @(x) x.*(-m.^2./q(x).^2 + (2.*m.*n)./q(x) - n.^2 + q(x).^(-2) - (2.*n)./(q(x).*m) + n.^2./m.^2 + (iar(x).*alpha(x).*(n.^4./m.^4 - n.^2./m.^2)) ) ; %%
    DRX_1D = @(x) x.*(3.*x.*(q(x).^(-2) - (2.*n)./(m.*q(x)) + n.^2./m.^2) - (2.*x.^2.*qp(x))./q(x).^3 + (2.*x.^2.*qp(x).*n)./(q(x).^2.*m) ) ;   %% le dernier + vient d etre changé avant c etait moins
    DRX_2D = @(x) x.*(x.^2.*(q(x).^(-2) - (2.*n)./(q(x).*m) + n.^2./m.^2) ) ;   %%%

    % Chi (= resistive correction) contribution to the main harmonic governing equation:    
    DRC_0D = @(x) x.*((3.*x.*n.*qp(x))./(m.*q(x).^2) + (x.^2.*n.*qpp(x))./(m.*q(x).^2) - (2.*x.^2.*n.*qp(x).^2)./(m.*q(x).^3)  - m.*n./q(x) + n.^2 + n./(m.*q(x)) - n.^2./m.^2)  ; %  vient
    DRC_1D = @(x) 3.*x.^2.*( n./(m.*q(x)) - n.^2./m.^2 ) - (alpha(x).*n.^2.*deltap(x).*x.^2)./(m.^2) ; %
    DRC_2D = @(x) x.*(x.^2.*( n./(m.*q(x)) - n.^2./m.^2 )) ; %%  % 
    
    % \Delta \Xi_{\Gamma} (=compressibility correction) contribution to the
    % main harmonic governing equation
    DRXG_0D = @(x) x.*((iar(x).*alpha(x).*n.^2./m.^2).*((n.^2./m.^2) - 1)) - x.*((n.*alpha(x).*deltap(x).*x.*qp(x))./(m.*q(x).^2)); 
    DRXG_1D = @(x) x.*(-1.*alpha(x).*n.^2.*deltap(x).*x.*(1 - m./(n.*q(x)))./m.^2) ;  %%
    DRXG_2D = @(x) 0.*x ; 
        
% -------------------------------------------------------------------------
    % Coefficients in Ohm's law (no difference wrt the paper or the
    % documentation):
% -------------------------------------------------------------------------
    
    % Xi contribution to Ohm's law:
    OLC_0D = @(x) x.^3 ;  %
    OLC_1D = @(x) 0.*x  ;
    OLC_2D = @(x) 0.*x  ;
        
    ORX_0D = @(x) (3.* Dqp(x).*x.^2 + x.^3.*Dqpp(x) + Dq(x).*(1-m.^2).*x ).*(rStar.^2/SL) ; %% 
    ORX_1D = @(x) (2.*Dqp(x).*x.^3 + 3.*Dq(x).*x.^2).*(rStar.^2/SL)   ; 
    ORX_2D = @(x) Dq(x).*x.^3.*(rStar.^2/SL) ;   
    
    % Chi contribution to Ohm's law:
    ORC_0D = @(x) (1-m.^2).*x.*(rStar.^2/SL)  ; 
    ORC_1D = @(x) 3.*x.^2.*(rStar.^2/SL)  ;   
    ORC_2D = @(x) x.^3.*(rStar.^2/SL)  ; 

% -------------------------------------------------------------------------
    % Coefficients in the definition of the compression variable 
    % \Delta \Xi_{\Gamma} (no difference wrt the paper or the documentation):
% -------------------------------------------------------------------------
    
    % \Delta \Xi_{\Gamma} contribution: 
    DXG_LXG_0D = @(x) (alfven(x).^2.*q(x).*(q(x) - m./n))./n.^2 ; % vient
    DXG_LXG_1D = @(x) 0.*x   ;
    DXG_LXG_2D = @(x) 0.*x  ;
        
    DXG_RXG_0D = @(x) (-1.*sound(x).^2.*(q(x)-m./n).^3)./q(x)  ;%%
    DXG_RXG_1D = @(x) 0.*x   ;
    DXG_RXG_2D = @(x) 0.*x   ;
    
    % Chi contribution:
    DXG_RC_0D = @(x) -1.*sound(x).^2.*(q(x)-m./n).^2  ;%%
    DXG_RC_1D = @(x) 0.*x   ;
    DXG_RC_2D = @(x) 0.*x   ;

end

function applyBCs()
    for ii=[1,N]
        DLX(ii,:) = 0; DLX(:,ii) = 0; DLX(ii,ii) = 0.;
        DRX(ii,:) = 0; DRX(:,ii) = 0; DRX(ii,ii) = 1.;
        DRC(ii,:) = 0; DRC(:,ii) = 0; DRC(ii,ii) = 0.;
        DRXG(ii,:) = 0; DRXG(:,ii) = 0; DRXG(ii,ii) = 0.;
        ORC(ii,:) = 0; ORC(:,ii) = 0; ORC(ii,ii) = 1.;
        OLC(ii,:) = 0; OLC(:,ii) = 0; OLC(ii,ii) = 0.;
        ORX(ii,:) = 0; ORX(:,ii) = 0; ORX(ii,ii) = 0.;
        DXG_LXG(ii, :) = 0; DXG_LXG(:, ii) = 0; DXG_LXG(ii, ii) = 0.;
        DXG_RXG(ii, :) = 0; DXG_RXG(:, ii) = 0; DXG_RXG(ii, ii) = 1.;
        DXG_RC(ii, :) = 0; DXG_RC(:, ii) = 0; DXG_RC(ii, ii) = 0.;
    end
    if m==1 %in this case Neumann BC's at r=0 for xi
        DRX(1, 2) = -1.;
        DLX(2, :) = 0; DRC(2,:) = 0; ORX(2,:) = 0; DRXG(2,:) =0;
        DRX(2,:) = -0; DRX(2,1) = -1; DRX(2,3) = 1.;
        ORC(1,2) = -1.;
        OLC(2, :) = 0; ORX(2,:) = 0.; DXG_RC(2,:) = 0.;
        ORC(2,:) = 0; ORC(2,1) = -1; ORC(2,3) = 1.;
        DXG_RXG(1,2) = -1.; 
        DXG_LXG(2,:) = 0; DXG_RC(2,:) = 0.; DRXG(2,:) = 0; 
        DXG_RXG(2,:) = 0.; DXG_RXG(2,1) = 1.; DXG_RXG(2,3) = -1.;
    end
end

function disc = discretize_S(fun, fund, fundd)
    ooo = ones(N-1,1)./dx ;
    op = [ooo;0] ; 
    om = [0;ooo] ; 
    Dm = ([-2;dx]+[dx;1])/2 ;
    D1F = sparse(diag(fund(xall)) + diag(fund(xall(1:end-1)),1) + diag(fund(xall(2:end)),-1)) ;
    D2F = sparse(diag(fundd(xall)) + diag(fundd(xall(1:end-1)),1) + diag(fundd(xall(2:end)),-1)) ;
    C0D = sparse(diag(fun(xall)));
    C1D = D1F.*sparse((diag(ooo,1) - diag(op)) + (diag(om) - diag(ooo,-1)))/2;
    C1D(1,1) = 2*C1D(1,1); C1D(1,2) = 2*C1D(1,2);
    C1D(end,end-1) = 2*C1D(end,end-1); C1D(end,end) = 2*C1D(end,end);
    C2D = D2F.*sparse(((diag(ooo./Dm(1:end-1),1) - diag(op./Dm)) - (diag(om./Dm) - diag(ooo./Dm(2:end),-1))));
    C2D(1,2) = C2D(1,2); C2D(end,end-1) = C2D(end,end-1);
    disc = C0D + C1D + C2D;
end

function makePlots()

    cl = EPFLcolors();  % EPFL color palette (in utils)
    colortab = {'groseille', 'rose', 'zinzolin', 'chartreuse', ...
        'montRose', 'ardoise', 'vertDeau'};
    colormat = cellfun(@(c) cl.(c), colortab, 'UniformOutput', false);
    colormat = cat(1, colormat{:});

    figure
    set(gca,'ColorOrder',colormat,'nextplot','replacechildren')
    set(0,'defaultTextInterpreter','tex');
    set(0, 'DefaultLineLineWidth', 3);
    set(groot,'defaultLineMarkerSize',12);
    set(gca,'TickLabelInterpreter','tex');
    set(0, 'DefaultAxesFontSize', 20)


    plot(xall, Xi, xall, Chi, xall, XiG, 'LineWidth',3,'MarkerSize',12);
    for i = 1:max(size(RS))
        xline(double(RS(i)), '--', 'LineWidth', 2);
    end
    ylabel('[a.u.]')
    xlabel('r / a')
    legend('\xi^r', [num2str(factorChi),' \chi'], [num2str(factorChi),' \Delta\xi_{\Gamma}'], 'r_s');
    xlim([0,1])
    titre = ['Resistive IM solver, \gamma / \omega_A =  ',num2str(ev),' \alpha = ', num2str(alpha(rStar)), ' s = ',' s = ', num2str(result.s)];
    title(titre,'interpreter','tex', 'FontWeight', 'normal')
    fontname('Times New Roman')
end

end