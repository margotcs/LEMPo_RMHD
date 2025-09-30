function [ev, result] = LEMPo_res(m, n, N, ev_guess, toPlot, profiles, opts)

% Run the full resistive model (so including infernal corrections)
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

resistiveSimu = 1;
input = getProfiles(m, n, resistiveSimu, profiles, toPlot, opts);
[q, qp, qpp, beta, betap, betapp, alpha, alphap, R0, B0, Gamma, rStar, RS, ...
    RSp, RSm, elong, triangByR] = deal(input.q, input.qp, input.qpp, input.beta, input.betap, input.betapp, ...
         input.alpha, input.alphap, input.R0, input.B0, input.Gamma, input.rStar, input.RS, input.RSp, ...
         input.RSm, input.elong, input.triangByR);

% Options and simulation flags
[additionalPlot,  addShaping, higherOscModes, reverse_qprofile,  ...
 NeumannBC, N_tries, SL, resOnlyOnMainRatSurf, neglectDeltaBParallel, ...
 noCompr] = deal(input.additionalPlot, input.addShaping, ...
         input.higherOscModes, input.reverse_qprofile, input.NeumannBC, ...
         input.N_tries, input.SL, input.resOnlyOnMainRatSurf, ...
         input.neglectDeltaBParallel, input.noCompr);

% Set up grid (bunched around rational surfaces)
xall = gridBunching(RS, RSp, RSm, N, input.upperBound, opts);
dx = diff(xall);

%--------------------------------------------------------------------------
%  Define coefficients for Chi variables
%--------------------------------------------------------------------------
% Since we use chi as the resistive variable for both for the main mode 
% and sideband, we need a coefficient to tell us when does it stand for the
% main harmonic contribution \chi^{(m)}, or the sideband contribution 
% chi^{(m+1)}. (See Appendix E of paper 'Fundamental properties of ideal 
% and resistive infernal modes in tokamaks', PPCF 2024, for more details).
% We will need these coefficients when discretizing the matrices.
% coefResMainChi will thus be 1 around the main rational surface. Around 
% the upper sideband rational surface: coefResMainChi=0, and coefResSidChi =1.


% In the current version MakeChiCoefNoSidRes, resistivity is actually included
% only for the main harmonic. coefResSidChi is always zero.
[coefResMainChi, coefResSidChi] = MakeChiCoefNoSidRes(m, xall,RS, input.upperBound);

% Need to pass the coefficients to the input struct for the makeFirstPlots
% function:
input.coefResMainChi = coefResMainChi;
input.coefResSidChi = coefResSidChi; 
%--------------------------------------------------------------------------


% The line below corresponds to higher order terms that we might want to
% turn on and off, to look at their impact. By default everything is taken
% into account (value 1).
AddS1 = 1;  AddS2 = 1; AddW1 = 1;  AddW2 = 1;  AddR1 = 1;  AddR2 = 1;
AddJ = 1;   
AddM = 1; AddR = 1; % these 2 appear only for m=1


% Multiplication factors for the variables on the plots
factorChi = 10;
factorSid = 1;
factorDeltaXi = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%      Definitions of Some Quantities and First Outputs 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

qs = m/n;
qsp = (m+1)/n;
qsm = (m-1)/n;
s = @(x) x./q(x).*qp(x);      % shear
iar = @(x) x./R0;             % inverse aspect ratio
% below is radial derivative of Shafranov shift
deltap = @(x) q(x).^2./((x+1e-10).^3).*cumtrapz(x,x.^3./(q(x).^2.*R0) + x.^2.*alpha(x)./(q(x).^2));  
dqqs = @(x) (q(x) - qs)./qs;
dqqsp = @(x) qp(x)./qs;

OneOverqp = @(x) -qp(x)./(q(x).^2); 
OneOverqpp = @(x) - qpp(x)./q(x).^2 + 2.*qp(x).^2./q(x).^3; 

Q = @(x) (q(x).^(-1) - 1/qs);
Qp = @(x) OneOverqp(x);
QPlus = @(x) (q(x).^(-1) - qsp^(-1));
QPlusp = @(x) OneOverqp(x);
QPluspp = @(x) OneOverqpp(x);
QMinus = @(x) (q(x).^(-1) - qsm^(-1));
QMinusp = @(x) OneOverqp(x);
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
result.RSp = RSp;
result.alphaRS = alpha(RS);
result.s = s(RS);
result.q = q(xall);
result.deltap = deltap(xall);
result.elong = elong;
result.triangByR = triangByR;
result.SL = SL;
result.omegaS = sound(xall);
result.omegaA = alfven(xall);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%      Plot Profiles (if requested)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if additionalPlot
    makeFirstPlots(input, xall, N, resistiveSimu) % (utils function)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%      Initialise & Define Coefficients in the Equations 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MG_0D = []; MG_1D = []; MG_2D = [];
M_m_0D = []; M_m_1D = []; M_m_2D = [];
M_sp_0D = []; M_sp_1D = []; M_sp_2D = [];
M_sm_0D = []; M_sm_1D = []; M_sm_2D = [];
M_chi_0D = []; M_chi_1D = []; M_chi_2D = [];
M_compr_0D = []; M_compr_1D = []; M_compr_2D = [];

SpG_0D = []; SpG_1D = []; SpG_2D = [];
Sp_sp_0D = []; Sp_sp_1D = []; Sp_sp_2D = [];
Sp_m_0D = []; Sp_m_1D = []; Sp_m_2D = [];
Sp_chi_0D = []; Sp_chi_1D = []; Sp_chi_2D = [];
Sp_chip_0D = []; Sp_chip_1D = []; Sp_chip_2D = [];
Sp_compr_0D = []; Sp_compr_1D = []; Sp_compr_2D = [];

SmG_0D = []; SmG_1D = []; SmG_2D = [];
Sm_sm_0D = []; Sm_sm_1D = []; Sm_sm_2D = [];
Sm_m_0D = []; Sm_m_1D = []; Sm_m_2D = [];
Sm_chi_0D = []; Sm_chi_1D = []; Sm_chi_2D = [];
Sm_compr_0D = []; Sm_compr_1D = []; Sm_compr_2D = [];

OG_0D = []; OG_1D = []; OG_2D = [];
O_m_0D = []; O_m_1D = []; O_m_2D = [];
O_chi_0D = []; O_chi_1D = []; O_chi_2D = [];
Op_sp_0D = []; Op_sp_1D = []; Op_sp_2D = [];
Op_chip_0D = []; Op_chip_1D = []; Op_chip_2D = [];

ComprG_0D = []; ComprG_1D = []; ComprG_2D = [];
Compr_compr_0D = []; Compr_compr_1D = []; Compr_compr_2D = [];
Compr_chi_0D = []; Compr_chi_1D = []; Compr_chi_2D = [];


defineLEMPoResOperators; % nested function which defines all operators 
%                         initialised above, model-specific.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%      Discretize Operators
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simple 1st order finite difference discretisation

MG = discretize_S(MG_0D,MG_1D,MG_2D);
M_m = discretize_S(M_m_0D,M_m_1D,M_m_2D);
M_chi = discretize_S(M_chi_0D,M_chi_1D,M_chi_2D).*coefResMainChi;
M_sp = discretize_S(M_sp_0D,M_sp_1D,M_sp_2D);
if m~=1
    M_sm = discretize_S(M_sm_0D, M_sm_1D, M_sm_2D);
end
if ~noCompr
    M_compr = discretize_S(M_compr_0D, M_compr_1D, M_compr_2D);
end

SpG = discretize_S(SpG_0D, SpG_1D, SpG_2D); 
Sp_sp = discretize_S(Sp_sp_0D, Sp_sp_1D, Sp_sp_2D);
Sp_m = discretize_S(Sp_m_0D, Sp_m_1D, Sp_m_2D);
if ~resOnlyOnMainRatSurf
    Sp_chi = discretize_S(Sp_chi_0D, Sp_chi_1D, Sp_chi_2D).*coefResMainChi + ...
        discretize_S(Sp_chip_0D, Sp_chip_1D, Sp_chip_2D).*coefResSidChi;
else % this case is the case where we dont have resistivity on the sideband
    Sp_chi = discretize_S(Sp_chi_0D, Sp_chi_1D, Sp_chi_2D).*coefResMainChi ;
end
if ~noCompr
    Sp_compr = discretize_S(Sp_compr_0D, Sp_compr_1D, Sp_compr_2D);
end

if m~=1
    SmG = discretize_S(SmG_0D, SmG_1D, SmG_2D);
    Sm_sm = discretize_S(Sm_sm_0D, Sm_sm_1D, Sm_sm_2D);
    Sm_m = discretize_S(Sm_m_0D, Sm_m_1D, Sm_m_2D);
    Sm_chi = discretize_S(Sm_chi_0D, Sm_chi_1D, Sm_chi_2D).*coefResMainChi;
    if ~noCompr
        Sm_compr = discretize_S(Sm_compr_0D, Sm_compr_1D, Sm_compr_2D);
    end
end

% Line below doesnt have a coefChi to set Chi to zero when no equation
% defines it.
OG = discretize_S(OG_0D,OG_1D,OG_2D);
O_m = discretize_S(O_m_0D,O_m_1D,O_m_2D).*coefResMainChi;
if ~resOnlyOnMainRatSurf
    O_chi = discretize_S(O_chi_0D,O_chi_1D,O_chi_2D).*coefResMainChi + ...
            discretize_S(Op_chip_0D,Op_chip_1D,Op_chip_2D).*coefResSidChi;
    O_sp = discretize_S(Op_sp_0D,Op_sp_1D,Op_sp_2D).*coefResSidChi;
else % i.e. if there is no resistivity on sidebands :
    O_chi = discretize_S(O_chi_0D,O_chi_1D,O_chi_2D).*coefResMainChi;
    O_sp = 0.*xall;
end

if ~noCompr
    ComprG = discretize_S(ComprG_0D, ComprG_1D, ComprG_2D);
    Compr_compr = discretize_S(Compr_compr_0D, Compr_compr_1D, Compr_compr_2D);
    Compr_chi = discretize_S(Compr_chi_0D, Compr_chi_1D, Compr_chi_2D).*coefResMainChi;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%      BCs (model-specific), and Build Matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nested functions 

applyBCs; 
[A,B] = buildLEMPoResMatrices();  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%      Solve Generalized Eigenvalue Problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[V, ev] = tryEigs(A, B, ev_guess, N_tries, xall, N, rStar, reverse_qprofile, ...
   higherOscModes, additionalPlot); % (utils function)

if isempty(V) % if couldnt get any growth rates: stop the execution
   ev = 0;
   return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%      Recover Eigenvalue and Normalise Eigenvectors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if sum(V(1:N)) < 0
    V(:) = -V(:);
end

Xi =  V(1:N) ;
maxXi = max(Xi);
Xi = Xi/maxXi;
Chi = V(2*N+1:3*N) ; Chi = factorChi.*Chi./maxXi;
Xip =  V(3*N+1:4*N) ; Xip = factorSid.*Xip./maxXi;

if m~=1
    if NeumannBC
        Xim =  V(4*N+1:5*N) ; Xim = factorSid.*Xim./maxXi;
        if ~noCompr
            DeltaXi = V(5*N+1: 6*N); DeltaXi = factorDeltaXi.*DeltaXi./maxXi;
        end
    else
         Xim =  V(5*N+1:6*N) ; Xim = factorSid.*Xim./maxXi;
         if ~noCompr
            DeltaXi = V(7*N+1: 8*N); DeltaXi = factorDeltaXi.*DeltaXi./maxXi;
         end
     end
else
    if ~noCompr
        if NeumannBC
            DeltaXi = V(4*N+1: 5*N); DeltaXi = factorDeltaXi.*DeltaXi./maxXi;
        else
            DeltaXi = V(5*N+1: 6*N); DeltaXi = factorDeltaXi.*DeltaXi./maxXi;
        end
    end
end

%% Go back to 'true' sidebands

 % We computed everything for Xip bar, so variable Xip is actually
 % the bar variable for now (Respectively Xim bar, and Xim for the lower
 % sideband). So the bit below calls the subfunction reverseChangeVariables
 % to return to the original sidebands expression. (See below eq. 10 in
 % paper "Fundamental properties of ideal and resistive infernal modes in
 % tokamaks")

 Xip = reverseChangeVariablesSid(xall, Xi, Chi, Xip,N,  "upper",ev,  m, dqqs,...
                deltap, rStar, RSp, NeumannBC, resistiveSimu, additionalPlot, ...
                coefResMainChi);
 if m~=1
     Xim = reverseChangeVariablesSid(xall, Xi, Chi, Xim,N,  "lower",ev, m, dqqs, ...
            deltap, rStar, RSp, NeumannBC, resistiveSimu, additionalPlot, ...
            coefResMainChi);
 end
% (utils function)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%      Build Result Structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


result.Xi = Xi;
result.ev = ev; % in the resistive code we obtain the growth rate non square,
% so no need for square root here
result.Xip = Xip;
result.Chi = Chi./factorChi;
if ~noCompr
    result.DeltaXi = DeltaXi./factorDeltaXi;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%      Plot Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if toPlot
    makePlots();
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%      Nested Functions (model-specific)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function defineLEMPoResOperators()
    % How are the elements named in this subfunction :
    %   X_Y_ZD : 
    %   - X corresponds to the equation we are looking at. It can be the main
    % harmonic equation ( 'M' ), the upper/lower sideband eq. ( 'Sp'/'Sm'),
    % or the Ohm's law if we add resistivity. 
    % Whenever the G letter appears, it means we are writing the term that
    % multiplies gamma^2 : it will be in the B matrix, when solving B  Xi  ev =
    % A Xi. All other terms will end up in the A matrix.
    %   - Y corresponds to a specific term in the equation, It can be term 
    % concerning the main harmonic displacement ('m'), or sidebands ('sp' or
    % 'sm')
    %   - Z is 1 if we are looking at the first derivative of the variable
    %   given by Y, 2 for second derivative, etc.
    
    % Example : write the main governing eq. like 
    % MG * (gamma^2/omega_A^2) = M_m_0D * xi^m + M_m_1D * xi^m ' + M_m_2D * 
    % xi^m '' + M_sp_0D * xi^m+1 + M_sp_1D * xi^m+1 ' + M_sp_2D * xi^m+1 ''
    % + M_sm_0D * xi^m-1 etc ...

    % Everything is multiplied by r in the code, wrt the paper or the
    % documentation.

% -------------------------------------------------------------------------
    % Coefficients in the main mode governing equation (they have an additional 
    % factor (-r) wrt the paper or the documentation):
% -------------------------------------------------------------------------
    
    % Main harmonic contribution to the main harmonic governing equation:
    MG_0D = @(x) -x.*(m^2-1).*(1+2*qs^2)./(m^2);
    MG_1D = @(x) 3.*x.^2.*(1+2*qs^2)./(m^2);
    MG_2D = @(x) x.^3.*(1+2*qs^2)./(m^2);
    
    if m==1
        M_m_0D = @(x) x.*(m^2-1).*Q(x).^2 - x.*alpha(x).*((x.*(qs^(-2) -1)./R0) - alpha(x)./4)./(qs^2);
        M_m_1D = @(x) -3.*(x.^2).*Q(x).^2 - (x.^3).*2.*Q(x).*Qp(x); 
        if AddM
            M_m_0D = @(x) M_m_0D(x) - x.*alpha(x).*iar(x).*dqqs(x)./(2.*qs^2);
        end
    else
        M_m_0D = @(x) x.*(m^2-1).*Q(x).^2 - alpha(x).*x.*(x.*(qs^(-2) -1)./R0 - alpha(x)./2)./(qs^2);%
        M_m_1D = @(x) -3.*(x.^2).*Q(x).^2 - (x.^3).*2.*Q(x).*Qp(x);
    end
    if neglectDeltaBParallel
        M_m_0D = @(x) M_m_0D(x) + x.*(alpha(x).^2)./(4*qs^4);
    end
    if addShaping %this is the same for all m 
        M_m_0D = @(x) M_m_0D(x) - x.*alpha(x).*((3/4).*(elong-1).*(1-2*triangByR*R0).*x./R0)./(qs^2);
    end
    if AddW1  
        if m==1
            factorW1 = @(x) dqqs(x).*qs.^(-2).*(4.*iar(x).^2.*(2 - qs.^(-2)) + iar(x).*alpha(x).*(5 - 2.*qs.^(-2)) ...
                            + (3/2).*alpha(x).^2 + deltap(x).*(12.*deltap(x) - (13/2).*alpha(x) ...
                            + (1/2).*x.*alphap(x) - 6.*iar(x)) ); %
            M_m_0D =  @(x) M_m_0D(x) - x.*factorW1(x);
        else
            factorW1 = @(x) dqqs(x).*qs.^(-2).*(4.*iar(x).^2.*(2 - qs.^(-2)) + iar(x).*alpha(x).*(5 - 2.*qs.^(-2)) ...
                            + 2.*alpha(x).^2 + deltap(x).*(12.*deltap(x) - 7.*alpha(x) + x.*alphap(x) - 6.*iar(x)) );
            M_m_0D =  @(x) M_m_0D(x) - x.*factorW1(x);%
        end
    end
    M_m_2D = @(x) -(x.^3).*Q(x).^2;
    
    % resistive contribution to the main mode equation :
    M_chi_0D = @(x) -3.*x.^2.*qp(x)./(qs.*q(x).^2) + x.^3.*(2.*qp(x).^2 ...
                    - qpp(x).*q(x))./(q(x).^3.*qs) + x.*(m^2 - 1).*Q(x)./qs;
    if AddJ
        if m==1
            factorJ = @(x) qs.^(-2).*(4.*iar(x).^2.*(2 - qs.^(-2)) + iar(x).*alpha(x).*((7/2) + 1.*qs.^(-2)) ...
                     + (3/2).*alpha(x).^2 + deltap(x).*(12.*deltap(x) - 7.*alpha(x) + x.*alphap(x) - 6.*iar(x)) ); %
        else
            factorJ = @(x) qs.^(-2).*(4.*iar(x).^2.*(2 - qs.^(-2)) + iar(x).*alpha(x).*(3 + 1.*qs.^(-2)) ...
                      + alpha(x).^2 + deltap(x).*(12.*deltap(x) - 7.*alpha(x) + x.*alphap(x) - 6.*iar(x)) );
        end
        M_chi_0D = @(x) M_chi_0D(x) + x.*factorJ(x);
    end
    M_chi_1D = @(x) -3.*x.^2.*Q(x)./qs; 
    if AddW2
        if m~=1
            M_m_1D = @(x) M_m_1D(x) - x.^2.*dqqs(x).*qs.^(-2).*deltap(x).*alpha(x);
            M_chi_1D = @(x) M_chi_1D(x) + x.^2.*deltap(x).*alpha(x)./qs.^2;
        end
    end
    if (m==1) && AddR
        M_chi_1D = @(x) M_chi_1D(x) + x.^2.*alpha(x).*deltap(x)./(qs.^2);
        M_m_1D = @(x) M_m_1D(x) - x.^2.*alpha(x).*deltap(x).*dqqs(x)./(2.*qs.^2);
    end
    M_chi_2D = @(x) -x.^3.*Q(x)./qs;
    
    % Sideband contribution to the main mode equation :
    M_sp_0D = @(x) -x.*alpha(x).*(2+m)./(2*(qs^2)*(1+m));
    M_sp_1D = @(x) -alpha(x).*(x.^2)./(2*(qs^2)*(1+m));
    M_sp_2D = @(x) 0.*x;
    if AddS1
       M_sp_0D = @(x) M_sp_0D(x) + (2+m).*x.*(dqqs(x).*(iar(x).*(1+m) + alpha(x) ...
                        + m.*deltap(x)))./((1+m).*qs.^2);
       M_sp_1D = @(x) M_sp_1D(x) + x.^2.*(dqqs(x).*( 2.*iar(x).*(1+m)   ...
                     + (2+m).*alpha(x) - (1+2.*m).*deltap(x)))./((1+m).*qs.^2);
    end
    if AddS2
        M_sp_2D = @(x) M_sp_2D(x) + x.^3.*dqqs(x).*deltap(x)./(qs.^2.*(1+m));
    end
    if m~=1 % lower sideband contribution to main mode eq.
       M_sm_0D = @(x) -x.*alpha(x).*(2-m)./(2.*(qs.^2).*(1-m));
       M_sm_1D = @(x) -alpha(x).*(x.^2)./(2.*(qs.^2).*(1-m));
       M_sm_2D = @(x) 0.*x; 
       if AddS1
           M_sm_0D = @(x) M_sm_0D(x) + (2-m).*x.*(dqqs(x).*(iar(x).*(1-m) + alpha(x) ...
                          - m.*deltap(x)))./((1-m).*qs.^2);
           M_sm_1D = @(x) M_sm_1D(x) + x.^2.*(dqqs(x).*( 2.*iar(x).*(1-m) ...
                         + (2-m).*alpha(x) - (1-2.*m).*deltap(x)))./((1-m).*qs.^2);
       end
       if AddS2     
           M_sm_2D = @(x) M_sm_2D(x) + x.^3.*dqqs(x).*deltap(x)./(qs.^2.*(1-m));
       end      
    end
    
    % Compression contribution to the equation

    if m==1  
       M_compr_0D = @(x) -x.*alpha(x).*((x.*(qs.^(-2) -1)./R0) - alpha(x)./4)./(qs^2);
       if AddM
           M_compr_0D = @(x) M_compr_0D(x) - x.*alpha(x).*iar(x).*dqqs(x)./(2.*qs.^2);
       end
    else
        M_compr_0D = @(x) (-x.*alpha(x).*(x.*(qs.^(-2) -1)./R0 - alpha(x)./2))./(qs.^(2));%
    end
    if AddJ  
        if m==1
            M_compr_0D = @(x) M_compr_0D(x) - x.*dqqs(x).*( 3.*iar(x).*alpha(x).*(0.5 - qs.^(-2)) ...
                            + 0.5.*deltap(x).*(alpha(x) - x.*alphap(x)))./qs.^2;
        else
            M_compr_0D = @(x) M_compr_0D(x) - x.*dqqs(x).*( iar(x).*alpha(x).*( 2- 3.*qs.^(-2)) ...
                             + alpha(x).^2)./(qs.^2);%
        end
    end
    if AddW2  
        if m==1
            M_compr_0D = @(x) M_compr_0D(x) + x.^2.*dqqsp(x).*deltap(x).*alpha(x)./(2.*qs.^2);%
        else
            M_compr_0D = @(x) M_compr_0D(x) + x.^2.*dqqsp(x).*deltap(x).*alpha(x)./(qs.^2);%
        end
    end
    M_compr_1D = @(x) 0.*x;  
    M_compr_2D = @(x) 0.*x;
    if m==1
        if AddR
            M_compr_0D = @(x) M_compr_0D(x) + x.^2.*alpha(x).*deltap(x).*dqqsp(x)./(2.*qs.^2);  %
            M_compr_1D = @(x) M_compr_1D(x) + x.^2.*alpha(x).*deltap(x).*dqqs(x)./(2.*qs.^2);   %
        end
    end


% -------------------------------------------------------------------------
    % Coefficients in the upper sideband governing equation
    % (they have an additional factor (-r) wrt the paper or the documentation):
% -------------------------------------------------------------------------
    
    % Upper sideband contribution to the upper sideband governing eq. :
    SpG_0D = @(x) -x.*m.*(m+2).*(1+2*qsp^2)./((m+1).^2);
    SpG_1D = @(x) 3.*x.^2.*(1+2.*qsp.^2)./((m+1).^2);
    SpG_2D = @(x) x.^3.*(1+2.*qsp.^2)./((m+1).^2);
    Sp_sp_0D = @(x) x.*m*(m+2).*QPlus(x).^2;
    Sp_sp_1D = @(x) -3.*(x.^2).*QPlus(x).^2 - (x.^3).*2.*QPlusp(x).*QPlus(x); 
    Sp_sp_2D = @(x) -(x.^3).*QPlus(x).^2;  %
    
    % Main mode and resistive contributions to the upper sideband
    % governing eq. :
    Sp_m_0D = @(x) x.*(alphap(x).*x - m.*alpha(x))./(2.*qs.^2.*(m+1));%
    Sp_m_1D = @(x) (x.^2).*alpha(x)./(2.*qs.^2.*(m+1));%
    Sp_m_2D = @(x) 0.*x; %
    Sp_chi_0D = @(x) 0.*x; % chi stands for the chi^{m} variable
    Sp_chi_1D = @(x) 0.*x;
    Sp_chi_2D = @(x) 0.*x;  %
        % chip stands for the chi^{m+1} variable
    Sp_chip_0D = @(x) 3.*x.^2.*OneOverqp(x)./qsp + x.^3.*OneOverqpp(x)./qsp ...
                      + x.*m.*(m+2).*QPlus(x)./qsp;
    Sp_chip_1D = @(x) -3.*x.^2.*QPlus(x)./qsp;
    Sp_chip_2D = @(x) -x.^3.*QPlus(x)./qsp;
    % compressibility correction to xip governing eq.
    Sp_compr_0D = @(x) x.*(alphap(x).*x - m.*alpha(x))./(2.*qs.^2.*(m+1)); %
    Sp_compr_1D = @(x) (x.^2).*alpha(x)./(2.*qs.^2.*(m+1));%
    Sp_compr_2D = @(x) 0.*x; %

    if AddR1
        Sp_m_0D = @(x) Sp_m_0D(x) + (-1).*(1+m).^(-1).*qs.^(-2).*x.*((1+m).*(((-1).*(3+m).*alpha(x) ...
                    +x.*alphap(x)+3.*(3+m).*deltap(x)).*dqqs(x)+x.*(alpha(x)+(-3).* ...
                    deltap(x)).*dqqsp(x))+((-2).*(1+m+m.^2).*dqqs(x)+(1+2.*m).*x.* ...
                    dqqsp(x)).*iar(x));
        Sp_m_1D = @(x) Sp_m_1D(x) + (-1).*(1+m).^(-1).*qs.^(-2).*x.^2.*dqqs(x).*((1+m).*(alpha(x)+( ...
                        -3).*deltap(x))+(1+2.*m).*iar(x));
        Sp_chi_0D = @(x) Sp_chi_0D(x) + (-1).*(1+m).^(-1).*qs.^(-2).*x.*((3+m.*(3+m)).*alpha(x)+(-1).* ...
                        m.*x.*alphap(x)+(-3).*(1+m).*(3+m).*deltap(x)+2.*(1+m+m.^2).*iar( ...
                        x)); 
        Sp_chi_1D = @(x) Sp_chi_1D(x) + (1+m).^(-1).*qs.^(-2).*x.^2.*(m.*alpha(x)+(-3).*(1+m).*deltap(x) ...
                            +iar(x)+2.*m.*iar(x));
        Sp_compr_0D = @(x) Sp_compr_0D(x) + (-1).*(1+m).^(-1).*qs.^(-2).*x.*(x.*alphap(x).*dqqs(x)+alpha(x) ...
                            .*((-1).*m.*dqqs(x)+x.*dqqsp(x)));               
        Sp_compr_1D = @(x) Sp_compr_1D(x) + (-1).*(1+m).^(-1).*qs.^(-2).*x.^2.*alpha(x).*dqqs(x);
    
    end 
    if AddR2
        Sp_m_0D = @(x) Sp_m_0D(x) - x.*( (2+m).*(iar(x) +alpha(x) -4.*deltap(x)).*dqqs(x))./(qs.^2); %
        Sp_chi_0D = @(x) Sp_chi_0D(x) + x.*((2+m).*(iar(x) +alpha(x) -4.*deltap(x)))./(qs.^2);   %
    end
 
% -------------------------------------------------------------------------
    % Coefficients in the lower sideband governing equation
    % (they have an additional factor (-r) wrt the paper or the documentation):
% -------------------------------------------------------------------------
 
    if m~=1 % if m=1 there is no lower sideband

        % Lower sideband contribution to the lower sideband governing eq.:
        SmG_0D = @(x) -x.*m.*(m-2).*(1+2.*qsm.^2)./((m-1)^2); %
        SmG_1D = @(x) 3.*x.^2.*(1+2.*qsm.^2)./((m-1)^2); %
        SmG_2D = @(x) x.^3.*(1+2.*qsm.^2)./((m-1)^2); %
    
        Sm_sm_0D = @(x) x.*m.*(m-2).*QMinus(x).^2;
        Sm_sm_1D = @(x) -3.*(x.^2).*QMinus(x).^2 - (x.^3).*2.*QMinusp(x).*QMinus(x);  
        Sm_sm_2D = @(x) -(x.^3).*QMinus(x).^2;
    
        % Main mode and resistive contributions to the lower sideband
        % governing eq. :
        Sm_m_0D = @(x) x.*(alphap(x).*x + m.*alpha(x))./(2.*qs.^2.*(1-m));
        Sm_m_1D = @(x) (x.^2).*alpha(x)./(2.*qs.^2.*(1-m));
        Sm_m_2D = @(x) 0.*x;
        Sm_chi_0D = @(x) 0.*x;
        Sm_chi_1D = @(x) 0.*x;
        Sm_chi_2D = @(x) 0.*x;
       
        Sm_compr_0D = @(x) x.*(alphap(x).*x + m.*alpha(x))./(2.*qs.^2.*(1-m));
        Sm_compr_1D = @(x) (x.^2).*alpha(x)./(2.*qs.^2.*(1-m));
        Sm_compr_2D = @(x) 0.*x;
        
        if AddR1
            Sm_m_0D = @(x) Sm_m_0D(x) + (-1).*qs.^(-2).*x.*((((-3)+m).*alpha(x)+x.*alphap(x)+(-3).*((-3) ...
                        +m).*deltap(x)).*dqqs(x)+x.*(alpha(x)+(-3).*deltap(x)).*dqqsp(x)+( ...
                        (-1)+m).^(-1).*(2.*(1+((-1)+m).*m).*dqqs(x)+((-1)+2.*m).*x.*dqqsp( ...
                        x)).*iar(x)); %
            Sm_m_1D = @(x) Sm_m_1D(x) + (-1).*((-1)+m).^(-1).*qs.^(-2).*x.^2.*dqqs(x).*(((-1)+m).*( ...
                            alpha(x)+(-3).*deltap(x))+((-1)+2.*m).*iar(x)); %            
            Sm_chi_0D = @(x) Sm_chi_0D(x) + ((-1)+m).^(-1).*qs.^(-2).*x.*((3+((-3)+m).*m).*alpha(x)+m.*x.* ...
                        alphap(x)+(-3).*((-3)+m).*((-1)+m).*deltap(x)+2.*(1+((-1)+m).*m).* ...
                        iar(x)); %              
            Sm_chi_1D = @(x) Sm_chi_1D(x) + (-1).*((-1)+m).^(-1).*qs.^(-2).*x.^2.*((-1).*m.*alpha(x)+3.*(( ...
                        -1)+m).*deltap(x)+iar(x)+(-2).*m.*iar(x)); %
                       
            Sm_compr_0D = @(x) Sm_compr_0D(x) + ((-1)+m).^(-1).*qs.^(-2).*x.*(x.*alphap(x).*dqqs(x)+alpha(x).*( ...
                            m.*dqqs(x)+x.*dqqsp(x)));
                    %
            Sm_compr_1D = @(x) Sm_compr_1D(x) + ((-1)+m).^(-1).*qs.^(-2).*x.^2.*alpha(x).*dqqs(x);
     
        end
        if AddR2
            Sm_m_0D = @(x) Sm_m_0D(x) - x.*( (2-m).*(iar(x) +alpha(x) -4.*deltap(x)).*dqqs(x))./(qs.^2); 
            Sm_chi_0D = @(x) Sm_chi_0D(x) + x.*( (2-m).*(iar(x) +alpha(x) -4.*deltap(x)))./(qs.^2); 
        end
    
    end

% -------------------------------------------------------------------------
    % Coefficients in Ohm's law (no difference wrt the paper or the
    % documentation):
% -------------------------------------------------------------------------

    OG_0D = @(x) x.^3 ;  % that is the coefficient of (gamma/omega_A) chi
    OG_1D = @(x) 0.*x  ;
    OG_2D = @(x) 0.*x  ;
    
    O_m_0D = @(x) (3.*Dqp(x).*x.^2 + x.^3.*Dqpp(x) + Dq(x).*(1-m.^2).*x ).*(rStar.^2./SL) ; %%
    O_m_1D = @(x) (2.*Dqp(x).*x.^3 + 3.*Dq(x).*x.^2).*(rStar.^2./SL)   ;
    O_m_2D = @(x) Dq(x).*x.^3.*(rStar.^2./SL) ;
    
    O_chi_0D = @(x) (1-m.^2).*x.*(rStar.^2./SL)  ;
    O_chi_1D = @(x) 3.*x.^2.*(rStar.^2./SL)  ;
    O_chi_2D = @(x) x.^3.*(rStar.^2./SL)  ;
      
    if ~resOnlyOnMainRatSurf % if resOnlyOnMainRatSurf, we would use Ohm law only
        % for the main mode \chi^{m}, so no need for the xi^{m+1} and chi^{m+1}
        % contributions that are respectively below 
    
        Op_sp_0D = @(x) (3.*QPlusp(x).*x.^2 + x.^3.*QPluspp(x) ...
                    + QPlus(x).*(1-(m+1).^2).*x ).*(qsp.*RSp.^2./SL) ; %
        Op_sp_1D = @(x) (2.*QPlusp(x).*x.^3 + 3.*QPlus(x).*x.^2).*(qsp.*RSp.^2./SL) ; %
        Op_sp_2D = @(x) QPlus(x).*x.^3.*(qsp.*RSp.^2./SL) ; %
        
        Op_chip_0D = @(x) (1-(m+1).^2).*x.*(qsp.*RSp.^2./SL) ;   %
        Op_chip_1D = @(x) 3.*x.^2.*(qsp.*RSp.^2./SL);  %;
        Op_chip_2D = @(x) x.^3.*(qsp.*RSp.^2./SL) ;  %
    end
    
% -------------------------------------------------------------------------
    % Coefficients in the definition of the compression variable
    % (no difference wrt the paper or the documentation):
% -------------------------------------------------------------------------

    ComprG_0D = @(x) (alfven(x).^2.*q(x).*(q(x) - m./n))./n.^2 ; 
    ComprG_1D = @(x) 0.*x;
    ComprG_2D = @(x) 0.*x;
        
    Compr_compr_0D = @(x) (-1.*sound(x).^2.*(q(x)-m./n).^3)./q(x)  ;
    Compr_compr_1D = @(x) 0.*x   ;
    Compr_compr_2D = @(x) 0.*x   ;
    
    Compr_chi_0D = @(x) -1.*sound(x).^2.*(q(x)-m./n).^2  ;
    Compr_chi_1D = @(x) 0.*x   ;
    Compr_chi_2D = @(x) 0.*x   ;


end

function applyBCs()

    for ii=[1,N]

        % B matrix, we solve B*V*ev = A*V
        MG(ii,:) = 0; MG(:,ii) = 0; MG(ii,ii) = 0.;
        SpG(ii, :)=0; SpG(:, ii) = 0; SpG(ii,ii) = 0.;
        if m~=1
            SmG(ii, :)=0; SmG(:, ii) = 0; SmG(ii,ii) = 0.;
        end
        OG(ii, :)=0; OG(:, ii) = 0; OG(ii,ii) = 0.;
        if ~noCompr
            ComprG(ii, :) = 0; ComprG(:,ii) = 0; ComprG(ii,ii) = 0.;
        end
        
        % A matrix: 
        % If m=1, need Neumann BC instead of Dirichlet on the main mode at
        % r=0
        M_m(ii,:) = 0; M_m(:,ii) = 0; 
        if m==1
             M_m(1,1) = 1.; M_m(1,2) = -1.;  M_m(N,N) = 1;
             M_m(2, :) = 0; O_m(2,:) = 0.; Sp_m(2,:) = 0.;
             MG(2,:) = 0.; M_chi(2,:) =0.; M_sp(2,:) =0.; 
             if ~noCompr
                M_compr(2,:) = 0.;
             end
             M_m(2,1) = -1.; M_m(2,3) = 1;
        else
            M_m(ii,ii) = 1.; 
        end
        M_chi(ii,:) = 0; M_chi(:, ii) = 0; M_chi(ii,ii) = 0;
        M_sp(ii,:) = 0; M_sp(:,ii) = 0; M_sp(ii,ii) = 0.;
        if ~noCompr
            M_compr(ii, :)= 0.; M_compr(:,ii) = 0.; M_compr(ii,ii) = 0.;
        end
        if m~=1
            M_sm(ii,:) = 0; M_sm(:,ii) = 0; M_sm(ii,ii) = 0.;
        end
        O_m(ii,:) = 0.; O_m(:,ii) = 0.; O_m(ii,ii) = 0.;
        O_sp(ii,:) = 0.; O_sp(:,ii) = 0.; O_sp(ii,ii) = 0.;

        %Neumann BC for chi, otherwise it doesnt behave correctly at r=0: 
        O_chi(ii, :) = 0; O_chi(:,ii) = 0; O_chi(ii, ii) =1;
        O_chi(1,2) = -1;
        O_sp(1,2) = -1; 
        O_chi(2, :) = 0; OG(2,:) = 0; Sp_chi(2,:) = 0;  M_chi(2,:)=0.;
        if ~noCompr
            Compr_chi(2,:) = 0;
        end
        O_m(2,:) = 0.;
        O_chi(2,1) = 1. ; O_chi(2,3) = -1; 

        O_sp(2,:) = 0; O_sp(2,1) = 1.; O_sp(2,3) = -1;

        % Neumann BC also for r=1, otherwise for some higher oscillating
        % modes the eigenvariables do not behave well.
        O_chi(N,N-1) = -1;
        O_sp(N,N-1) = -1;
        O_chi(N-1, ii) = 0; O_chi(N-1,N) = 1. ; O_chi(N-1,N-2) = -1;
        O_m(N-1,:) = 0; OG(N-1,:) = 0;
        Sp_chi(N-1,:) = 0;
 
        if NeumannBC
            Sp_sp(ii,:) = 0 ; Sp_sp(:,ii) = 0 ; Sp_sp(ii,ii) = 1.;
            Sp_sp(N,N) = 1.;
            Sp_sp(N, N-1) = -1; 
            Sp_sp(N-1, :) = 0.; Sp_sp(N-1, N) = 1.; Sp_sp(N-1, N-2) = -1; 
            Sp_chi(N-1,:) = 0.;  Sp_m(N-1,:) = 0.;
            O_chi(N,N-1) = -1;
            O_chi(N-1, ii) = 0; O_chi(N-1,N) = 1. ; O_chi(N-1,N-2) = -1;
            O_m(N-1,:) = 0; OG(N-1,:) = 0;
            
            M_sp(N-1,:) = 0.;
            M_chi(N-1,:) = 0.;
            if ~noCompr
                Sp_compr(N-1,:) = 0.;
            end
        else
            Sp_sp(ii,:) = 0.; Sp_sp(:,ii) = 0.; Sp_sp(ii,ii) = 1.;
        end
        Sp_m(ii,:) = 0 ; Sp_m(:,ii) = 0 ; Sp_m(ii,ii) = 0.;
        Sp_chi(ii,:) = 0; Sp_chi(:,ii) = 0; Sp_chi(ii,ii) = 0.;
        if ~noCompr
             Sp_compr(ii,:) = 0.; Sp_compr(:,ii) = 0.; Sp_compr(ii,ii) = 0.;
        end
        if m~=1
            Sm_sm(ii,:) = 0 ; Sm_sm(:,ii) = 0 ; 
            if m == 2   
                % Lower sid will be m = 1 : apply Neumann BC instead of 
                % Dirichlet 
                Sm_sm(1,1) = 1.; Sm_sm(1,2) = -1.;  
                Sm_sm(2,ii) = 0; Sm_sm(2,1) = 1.; Sm_sm(2,3) = -1.; Sm_sm(N,N) = 1;
            else
                Sm_sm(ii,ii) = 1.;
            end
            Sm_m(ii,:) = 0 ; Sm_m(:,ii) = 0; Sm_m(ii,ii) = 0. ; 
            Sm_chi(ii,:) = 0; Sm_chi(:,ii) = 0; Sm_chi(ii,ii) = 0.;
            if ~noCompr
                Sm_compr(ii,:) = 0; Sm_compr(:,ii) = 0; Sm_compr(ii,ii) = 0;
            end
        end
        if ~noCompr
            Compr_compr(ii,:) = 0; Compr_compr(:,ii) = 0; Compr_compr(ii,ii) = 1.;
            Compr_chi(ii,:) = 0; Compr_chi(:, ii) = 0; Compr_chi(ii,ii) =0.;
        end
    end

end

function [A, B] = buildLEMPoResMatrices()
    % This function builds A and B to solve A v = gamma/omegaA B v
    Z = sparse(N,N) ;
    S = speye(N,N);
    
    if NeumannBC % in this case, no sideband inertia, one less variable: gamma*xip
        SpG = Z;  
        if m==1 % if m=1, no lower sideband
          B = [ S, Z, Z,  Z,  Z,  Z
                Z,MG, Z,  Z,  Z,  Z
                Z, Z, OG, Z,  Z,  Z
                Z, Z, Z,  Z,  Z,  Z
                Z, Z, Z,  Z,  S,  Z
                Z, Z, Z,  Z,  Z, ComprG];

          A = [  Z,    S,   Z,     Z,           Z,        Z;
                M_m,   Z, M_chi, M_sp,         M_compr,   Z;
                O_m,   Z, O_chi,   Z,           Z,        Z;
               Sp_m,   Z, Sp_chi, Sp_sp,      Sp_compr,   Z; 
                 Z,    Z,    Z,     Z,          Z,        S;
                 Z,    Z,Compr_chi, Z,      Compr_compr,  Z]; 

        else
          B = [ S, Z, Z,  Z, Z, Z, Z,  Z
                Z,MG, Z,  Z, Z, Z, Z,  Z
                Z, Z, OG, Z, Z, Z, Z,  Z
                Z, Z, Z,  Z, Z, Z, Z,  Z
                Z, Z, Z,  Z, S, Z, Z,  Z 
                Z, Z, Z,  Z, Z,SmG,Z,  Z
                Z, Z, Z,  Z, Z, Z, S,  Z
                Z, Z, Z,  Z, Z, Z, Z, ComprG];
          
          A = [  Z,    S,   Z,     Z,     Z,   Z,        Z,        Z;
                M_m,   Z, M_chi, M_sp,  M_sm,  Z,       M_compr,   Z;
                O_m,   Z, O_chi,   Z,     Z,   Z,        Z,        Z;
               Sp_m,   Z, Sp_chi, Sp_sp,  Z,   Z,      Sp_compr,   Z; 
                 Z,    Z,   Z,     Z,     Z,   S,         Z,       Z;
               Sm_m,   Z, Sm_chi,  Z,    Sm_sm,Z,      Sm_compr,   Z;
                 Z,    Z,   Z,     Z,     Z,   Z,         Z,       S;
                 Z,    Z,Compr_chi,Z,     Z,   Z,    Compr_compr,  Z]; 
        end
    end
    
    if ~NeumannBC
        if m==1
            if noCompr % no need for the compression variable
                B = [ S, Z, Z, Z,  Z;
                      Z,MG, Z, Z,  Z;
                      Z, Z, OG, Z, Z;
                      Z, Z, Z,  S, Z;
                      Z, Z, Z, Z, SpG];
        
                A = [  Z, S,    Z,    Z,       Z;
                     M_m, Z, M_chi, M_sp,      Z;
                     O_m, Z, O_chi,   O_sp,    Z;
                      Z,  Z,    Z,     Z,      S;
                    Sp_m, Z, Sp_chi,  Sp_sp,   Z];
            else     
                B = [ S, Z, Z, Z,  Z,  Z,  Z
                      Z,MG, Z, Z,  Z,  Z,  Z
                      Z, Z, OG, Z, Z,  Z,  Z
                      Z, Z, Z,  S, Z,  Z,  Z
                      Z, Z, Z, Z, SpG, Z,  Z
                      Z, Z, Z, Z,  Z,  S,  Z
                      Z, Z, Z, Z,  Z,  Z, ComprG];
        
                A = [  Z, S,    Z,    Z,       Z,        Z,      Z;
                     M_m, Z, M_chi, M_sp,      Z,      M_compr,  Z;
                     O_m, Z, O_chi,   O_sp,    Z,        Z,      Z;
                      Z,  Z,    Z,     Z,      S,        Z,        Z;
                    Sp_m, Z, Sp_chi,  Sp_sp,   Z,      Sp_compr,   Z;
                      Z,  Z,    Z,     Z,      Z,       Z,          S;
                      Z,  Z,Compr_chi, Z,      Z,  Compr_compr,    Z];
            end
      else 
           if noCompr

           B = [ S, Z, Z, Z,  Z, Z,  Z;
                 Z,MG, Z, Z,  Z, Z,  Z; 
                 Z, Z, OG, Z, Z, Z,  Z; 
                 Z, Z, Z,  S, Z, Z,  Z; 
                 Z, Z, Z, Z, SpG,Z,  Z; 
                 Z, Z, Z, Z, Z,  S,  Z; 
                 Z, Z, Z, Z,  Z, Z, SmG];  

            
    
            A = [  Z,  S,     Z,    Z,    Z,    Z,    Z;
                  M_m, Z, M_chi,  M_sp,   Z, M_sm,    Z;   
                  O_m, Z, O_chi,  O_sp,   Z,    Z,    Z;    
                   Z,  Z,    Z,     Z,    S,    Z,    Z;    
                 Sp_m, Z, Sp_chi, Sp_sp,  Z,    Z,    Z;  
                   Z,   Z,    Z,    Z,    Z,    Z,    S;   
                 Sm_m, Z,  Sm_chi,  Z,    Z, Sm_sm,   Z]; 


           else

           B = [ S, Z, Z, Z,  Z, Z,  Z,  Z,  Z
                 Z,MG, Z, Z,  Z, Z,  Z,  Z,  Z
                 Z, Z, OG, Z, Z, Z,  Z,  Z,  Z
                 Z, Z, Z,  S, Z, Z,  Z,  Z,  Z
                 Z, Z, Z, Z, SpG,Z,  Z,  Z,  Z
                 Z, Z, Z, Z, Z,  S,  Z,  Z,  Z
                 Z, Z, Z, Z,  Z, Z, SmG, Z,  Z 
                 Z, Z, Z, Z,  Z, Z,  Z,  S,  Z
                 Z, Z, Z, Z,  Z, Z,  Z,  Z, ComprG];
            
    
            A = [  Z,  S,     Z,    Z,    Z,    Z,    Z,      Z,       Z;
                  M_m, Z, M_chi,  M_sp,   Z, M_sm,    Z,   M_compr,    Z;
                  O_m, Z, O_chi,  O_sp,   Z,    Z,    Z,      Z,       Z;
                   Z,  Z,    Z,     Z,    S,    Z,    Z,      Z,       Z;
                 Sp_m, Z, Sp_chi, Sp_sp,  Z,    Z,    Z,  Sp_compr,    Z;
                   Z,   Z,    Z,    Z,    Z,    Z,    S,      Z,       Z;
                 Sm_m, Z,  Sm_chi,  Z,    Z, Sm_sm,   Z, Sm_compr,     Z;
                   Z,  Z,     Z,    Z,  Z,    Z,     Z,     Z,         S;
                   Z,  Z,Compr_chi, Z,  Z,    Z,     Z,  Compr_compr, Z];
           end
        end
    end

    % To debug: check if matrices contain NaN values
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

end

function disc=discretize_S(fun,fund,fundd) % not model-specific but needs too many arguments
    % Simple 1D finite difference discretisation
    ooo = ones(N-1,1)./dx ;
    op = [ooo;0] ; 
    om = [0;ooo] ; 
    Dm = ([-2;dx]+[dx;1])/2 ;
    
    D1F = sparse(diag(fund(xall)) + diag(fund(xall(1:end-1)),1) + diag(fund(xall(2:end)),-1)) ;
    D2F = sparse(diag(fundd(xall)) + diag(fundd(xall(1:end-1)),1) + diag(fundd(xall(2:end)),-1)) ;
     
    C0D = sparse( diag( fun(xall) )) ;
    
    C1D = D1F.*sparse( (diag(ooo,1) - diag(op)) + (diag(om) - diag(ooo,-1)) )/2 ;
    C1D(1,1) = 2*C1D(1,1) ; C1D(1,2) = 2*C1D(1,2) ;
    C1D(end,end-1) = 2*C1D(end,end-1) ; C1D(end,end) = 2*C1D(end,end) ;
    
    C2D = D2F.*sparse( ((diag(ooo./Dm(1:end-1),1) - diag(op./Dm)) - (diag(om./Dm) - diag(ooo./Dm(2:end),-1))) ) ;
     C2D(1,2) = C2D(1,2); C2D(end,end-1) = C2D(end,end-1);
     
    
    disc = C0D + C1D + C2D ;
    
end

function makePlots()

    cl = EPFLcolors();  % EPFL color palette (in utils)
    colortab = {'groseille', 'canard', 'rose', 'zinzolin', 'chartreuse', ...
        'montRose', 'ardoise', 'vertDeau'};
    colormat = cellfun(@(c) cl.(c), colortab, 'UniformOutput', false);
    colormat = cat(1, colormat{:});

    figure
    set(0,'defaultTextInterpreter','tex');   
    set(0, 'DefaultLineLineWidth', 3);
    set(groot,'defaultLineMarkerSize',12);
    set(gca,'TickLabelInterpreter','tex');
    set(0, 'DefaultAxesFontSize', 20)
    set(gca,'ColorOrder',colormat,'nextplot','replacechildren')
 
    if m~=1
       plot(xall, Xi, '.-',xall, Xip, '.-', xall, Xim, '.-', xall, Chi, '.-','LineWidth',3,'MarkerSize',12);
       hold on
       if ~noCompr
            plot(xall, DeltaXi, '.-', 'LineWidth',3,'MarkerSize',12);
       end
       title(['m = ',num2str(m),' ,n= ',num2str(n),' ,\alpha(', num2str(rStar),...
             ') = ',num2str(alpha(rStar)),', \gamma / \omega_A =  ',num2str(ev)],...
           'interpreter','tex', 'FontWeight', 'normal');
       for i=1:max(size(RS))
           xline(double(RS(i)), '--', 'LineWidth', 2, 'Color','#B51F1F');
       end
       if ~isempty(RSp)
            xline(RSp, '--', 'LineWidth', 2,'Color', '#007480');
       end
       if ~noCompr
            legend('\xi^{(m)}', [num2str(factorSid),' \xi^{(m+1)}'], [num2str(factorSid),...
                 ' \xi^{(m-1)}'], [num2str(factorChi),' \chi'], [num2str(factorDeltaXi), ...
                    ' \Delta \xi_{\Gamma}'], 'r_s', 'r_{m+1}')
       else
           legend('\xi^{(m)}', [num2str(factorSid),' \xi^{(m+1)}'], [num2str(factorSid),...
                 ' \xi^{(m-1)}'], [num2str(factorChi),' \chi'], 'r_s', 'r_{m+1}')
       end
        
       ylabel('[a.u]');
       xlabel('r / a');       
    else
       plot(xall, Xi, '.-',xall, Xip, '.-', xall, Chi, '.-', 'LineWidth',3,'MarkerSize',12) ;
       hold on
       if ~noCompr
            plot(xall, DeltaXi, '.-', 'LineWidth',3,'MarkerSize',12);
       end
       title(['m = ',num2str(m),' ,n= ',num2str(n),' ,\alpha(', num2str(rStar),') = ',...
           num2str(alpha(rStar)), ', \gamma / \omega_A =  ',num2str(ev),', s(', num2str(rStar),') = ', ...
           num2str(s(rStar))],'interpreter','tex', 'FontWeight', 'normal');
        for i=1:max(size(RS))
           xline(double(RS(i)), '--', 'LineWidth', 2, 'Color','#B51F1F');
       end
       if ~isempty(RSp)
            xline(RSp, '--', 'LineWidth', 2,'Color', '#007480');
       end
       if ~noCompr
           legend('\xi^{(m)}', [num2str(factorSid),' \xi^{(m+1)}'], [num2str(factorChi),...
                  ' \chi'], [num2str(factorDeltaXi), ' \Delta \xi_{\Gamma}'], 'r_s', 'r_{m+1}')
       else
           legend('\xi^{(m)}', [num2str(factorSid),' \xi^{(m+1)}'], [num2str(factorChi),...
              ' \chi'], 'r_s', 'r_{m+1}')
       end
       ylabel('[a.u]');
       xlabel('r / a');
    end
    fontname('Times New Roman')
end


end