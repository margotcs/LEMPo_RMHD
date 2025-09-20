function [ev, result]=IM_LEMPo_ideal(m, n, N, ev_guess, toPlot, profiles, opts)

% Run the interchange ideal model (so without infernal corrections)
%
% INPUTS:
%   m, n                - Poloidal & toroidal mdoe numbers of the main
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

resistiveSimu = 0;
% Load profiles and physical parameters (assumes getProfiles is available)
input = getProfiles(m, n, resistiveSimu, profiles, toPlot, opts);
[q, qp, beta, betap, alpha, alphap, R0, rStar, RS] = deal( ...
    input.q, input.qp, input.beta, input.betap, ...
    input.alpha, input.alphap, input.R0, input.rStar, input.RS);
% Options:
[additionalPlot, addShaping, higherOscModes, reverse_qprofile, N_tries] = ...
    deal(input.additionalPlot, input.addShaping, input.higherOscModes, ...
    input.reverse_qprofile,  input.N_tries);


if addShaping
    error(['You set opts.addShaping but shaping of flux surfaces is not ',...
        'included in the Interchange Models.'])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%      Set up grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set up grid, bunched around the rational surface of the main harmonic
xall = gridBunching(RS, [], [], N, 1, opts); % (utils function)
dx = diff(xall);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%      Definitions of Some Quantities and First Outputs 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% radial derivative of Shafranov shift
deltap = @(x) q(x).^2./((x+1e-10).^3).*cumtrapz(x, x.^3./(q(x).^2.*R0) + x.^2.*alpha(x)./(q(x).^2));
s = @(x) x./q(x).*qp(x);    % shear
eps = @(x) x/R0;            % inverse aspect ratio

result.rgrid = xall;
result.alpha = alpha(xall);
result.alphaRS = alpha(RS);
result.s = s(rStar);
result.q = q(xall);
result.RS = RS;
result.mercier = alpha(RS).*RS.*((n/m)^2 -1)./R0 - s(RS)./4;
result.deltap = deltap(xall);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%      Plot Profiles (if requested)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if additionalPlot
    makeFirstPlots(input, xall, N, resistiveSimu); % (utils function)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%      Initialise & Define Coefficients in the Equations 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DLX_0D = []; DLX_1D = []; DLX_2D = [];
DRX_0D = []; DRX_1D = []; DRX_2D = [];

defineIMidealOperators; % nested function which defines operators 
%                         initialised above, model-specific

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%      Discretize Operators
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Simple 1st order finite difference discretisation
DLX = discretize_S(DLX_0D, DLX_1D, DLX_2D);
DRX = discretize_S(DRX_0D, DRX_1D, DRX_2D);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%      Boundary Conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
applyBCs; % nested function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%      Solve Eigenvalue Problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A = DRX; B = DLX;
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


[V, ev] = tryEigs(A, B, ev_guess, N_tries, xall, N, rStar, reverse_qprofile, ...
    higherOscModes, additionalPlot); % (utils function)

if isempty(V) % if couldnt get any growth rates: stop the execution
    ev = 0;
    return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%      Recover Eigenvalue, Normalise Eigenvectors, Build Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ev = sqrt(ev); % because for the ideal problem, what we obtain initially
% is gamma/omegaA squared.

if sum(V(1:N)) < 0
     V(:) = - V(:);
end 
Xi =  V(1:N) ;
maxXi = max(Xi);
Xi = Xi./maxXi; % radial plasma displacement of the main harmonic

result.ev = ev ;
result.Xi = Xi;
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

function defineIMidealOperators()
    % How are the elements named in this subfunction :
    %   X_ZD : 
    %   - X corresponds to the side of the dispersion relation 
    % we are looking at. DLX marks the terms which will multiply gamma^2
    % (inertia), it will be in the B matrix, when solving B  Xi  ev =
    % A Xi. DRX will be in the A matrix (contains Mercier term, FLB etc)
    %   - Z is 1 if we are looking at the first derivative of the variable
    %   given by Y, 2 for second derivative, etc.
    
    % Everything is multiplied by r in the code, wrt the paper or the
    % documentation.
    
    DLX_0D = @(x) x - x./m.^2 + (2.*x.*m.^2)./n.^2 - 2.*x/n.^2  ; %%%
    DLX_1D = @(x) -3.*x.^2.*m.^(-2) - 2.*x.^2./n.^2  ;    % 
    DLX_2D = @(x) -x.^3.*m.^(-2) - 2.*x.^3./n.^2  ;    %  
           
    DRX_0D = @(x) x.*(-m.^2./q(x).^2 + (2.*m.*n)./q(x) - n.^2 + q(x).^(-2) - (2.*n)./(q(x).*m) + n.^2./m.^2 + (eps(x).*alpha(x).*(n.^4./m.^4 - n.^2./m.^2)) ) ; %%
    DRX_1D = @(x) x.*(3.*x.*(q(x).^(-2) - (2.*n)./(m.*q(x)) + n.^2./m.^2) - (2.*x.^2.*qp(x))./q(x).^3 + (2.*x.^2.*qp(x).*n)./(q(x).^2.*m) ) ;   %% le dernier + vient d etre changé avant c etait moins
    DRX_2D = @(x) x.*(x.^2.*(q(x).^(-2) - (2.*n)./(q(x).*m) + n.^2./m.^2) ) ;   %%%
       
end

function applyBCs()
    for ii=[1,N]
        DLX(ii,:) = 0; DLX(:,ii) = 0; DLX(ii,ii) = 0.;
        DRX(ii,:) = 0; DRX(:,ii) = 0; DRX(ii,ii) = 1.;
    end
    if m==1 %in this case Neumann BCs at r=0 
        DRX(1, 2) = -1.;
        DLX(2, :) = 0; DRX(2,:) = -0; DRX(2,1) = -1; DRX(2,3) = 1.;
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
    figure
    set(0,'defaultTextInterpreter','tex');   
    set(0, 'DefaultLineLineWidth', 3);
    set(groot,'defaultLineMarkerSize',12);
    set(gca,'TickLabelInterpreter','tex');
    set(0, 'DefaultAxesFontSize', 20)

    plot(xall, Xi, '.-','LineWidth',3,'MarkerSize',12)
    for i = 1:max(size(RS))
       xline(double(RS(i)), '--b', 'LineWidth', 2);
    end
    xlim([0,1])
    ylabel('Radial plasma displacement [a.u.]')
    xlabel('r / a')
    legend('\xi^{(m)}', 'r_s');
    titre = ['Ideal IM solver, \gamma / \omega_A =  ',num2str(ev),' \alpha = ', num2str(alpha(rStar)), ' s = ', num2str(result.s)];
    title(titre,'interpreter','tex', 'FontWeight', 'normal','FontSize',15)
end

end