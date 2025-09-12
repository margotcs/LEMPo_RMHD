function input = getProfiles(m, n, resistiveSimu, profiles, toPlot, opts)
% Set up q, qp, qpp, beta, betap, betapp, alpha, alphap, R0, etc.
% Also updates all options using default values or user supplied through
% struct opts. Provides the complete input to run the model.

% INPUTS:
%   m, n                - Poloidal & toroidal mdoe numbers of the main
%                         harmonic you want to look at
%   resistiveSimu       - Boolean: 1 if this is a resistive model running
%   profiles            - Struct containing all necessary inputs, see
%                         documentation
%   toPlot              - Boolean: 1 if you want outputs plots, 0 otherwise
%   opts                - Struct containing options, can be set to []
%
% OUTPUTS:
%   input               - Complete input to run the model (struct)


input = profiles; % start with profiles, will add missing fields either 
%                   from default values or from opts. supplied.
input.m = m;
input.n = n;
x = sym('x');


% --- Physical constants 
if ~isempty(opts) && isfield(opts, 'B0') 
    input.B0 = opts.B0;
else
    input.B0=3; %default B0 (useful for resistive models) [T]
end
if ~isempty(opts) && isfield(opts, 'R0') 
    input.R0 = opts.R0;
else
    input.R0=3; %default major radius [m]
end
input.Gamma = 5/3; % Adiabatic index


% --- q profile and derivatives
qs = m/n;
if ~isempty(profiles) && isfield(profiles, 'q')
    input.q = profiles.q;
    if isfield(profiles, 'qp')
         input.qp = profiles.qp;
         if isfield(profiles, 'qpp')
            input.qpp = profiles.qpp;
         else
             error(['If you give an expression for derivative of q, you must', ...
                 'also give one for the second derivative : .qpp field.']);
         end
    else
        input.qp = eval(['@(x)' strrep(strrep(strrep(char(diff(input.q(x))),'^','.^'),'*','.*'),'/','./')]);
        input.qpp = eval(['@(x)' strrep(strrep(strrep(char(diff(input.qp(x))),'^','.^'),'*','.*'),'/','./')]);
    end
else % default q profile: resonant ultra flat.
    c = 32.33;
    d = 2.507;
    delta_q = 0.01;
    input.q = @(x) (qs - delta_q).*(1+ c.*x.^(2*d)).^(1/d);
    input.qp = eval(['@(x)' strrep(strrep(strrep(char(diff(input.q(x))),'^','.^'),'*','.*'),'/','./')]);
    input.qpp = eval(['@(x)' strrep(strrep(strrep(char(diff(input.qp(x))),'^','.^'),'*','.*'),'/','./')]);
end

% --- Find rational surface(s), define characteristic length rStar

% In the resistive solver we sometimes need a characteristic length
% (compute Lundquist number, etc), this is rStar. It will be the first
% rational surface encountered.

if isfield(profiles, 'RS')
   input.RS = profiles.RS;
   input.rStar = profiles.RS;
else
   input.RS = 0;
   for n_vpaSolve = 1:20
       temp = vpasolve(input.q(x) == m/n,x,[(n_vpaSolve-1)/20,n_vpaSolve/20],'Random', true);
       if ~isempty(temp)
            input.RS = [input.RS, temp];
       end
   end
   input.RS(input.RS == 0) = [];
   input.RS = unique(input.RS, 'sorted')
   if isempty(input.RS)
      error(['No rational surface could be found. Try changing the', ...
            'q profile, or supply a value for RS through profiles.RS']);
   end
   input.RS = double(input.RS)';
   input.rStar = double(input.RS(1)); 
end
  
% rational surfaces of the sidebands :
if isfield(profiles, 'RSp')
   input.RSp = profiles.RSp;
else
    try 
        RSp = vpasolve(q(x) == (m+1)/n , x ,[0,1]); 
        input.RSp = double(RSp);
        if isempty(nonzeros(input.RSp))
            input.RSp = [];
        end
    catch
        input.RSp = [];
    end
end


if isfield(profiles, 'RSm')
   input.RSm = profiles.RSm;
else
    try 
        RSm = vpasolve(q(x) == (m-1)/n , x ,[0,1]); 
        input.RSm = double(RSm);
        if isempty(nonzeros(input.RSm))
            input.RSm = [];
        end
    catch
        input.RSm = [];
    end
end

% --- Pressure profile, alpha, derivatives

if ~isempty(profiles) && isfield(profiles, 'beta')
    if isfield(profiles,'alpha') 
      error(['You want to specify both the beta profile and force', ...
     'the value of alpha at a certain radial location. Choose between the two.'])
    end
    input.beta = profiles.beta;
    if isfield(profiles, 'betap')
        input.betap = profiles.betap;
        input.betapp = profiles.betapp;
    else
        input.betap = eval(['@(x)' strrep(strrep(strrep(char(diff(input.beta(x))),'^','.^'),'*','.*'),'/','./')]);
        input.betapp = eval(['@(x)' strrep(strrep(strrep(char(diff(input.betap(x))),'^','.^'),'*','.*'),'/','./')]);
    end
    input.alpha = @(x) -(input.q(x).^2).*input.R0.*input.betap(x);
    input.alphap = @(x) -2.*input.qp(x).*input.q(x).*input.R0.*input.betap(x) -(input.q(x).^2).*input.R0.*input.betapp(x); 
else
    if ~isempty(profiles) && isfield(profiles, 'alpha') 
        alpha_wanted = profiles.alpha;
        if isfield(profiles, 'rSetAlpha') % value at which we set alpha
           rSetAlpha = profiles.rSetAlpha;
        end
    else
        alpha_wanted = 0.1;
    end
    if ~exist('rSetAlpha', 'var')
           rSetAlpha = input.rStar; 
    end
   
    if alpha_wanted >= 0 % set up a parabolic pressure profile
        
        beta0 = alpha_wanted/(2*input.q(rSetAlpha)^2*input.R0*rSetAlpha);
        input.beta = @(x) beta0*(1-x.^2).*(1-x.^15);
        input.betap = eval(['@(x)' strrep(strrep(strrep(char(diff(input.beta(x))),'^','.^'),'*','.*'),'/','./')]);
        input.betapp = eval(['@(x)' strrep(strrep(strrep(char(diff(input.betap(x))),'^','.^'),'*','.*'),'/','./')]);
        input.alpha = @(x) -(input.q(x).^2).*input.R0.*input.betap(x);
        input.alphap = @(x) -2.*input.qp(x).*input.q(x).*input.R0.*input.betap(x) -(input.q(x).^2).*input.R0.*input.betapp(x); 

    else % hollow pressure profile case
        r0 = 0.4; % r from which the pressure profile decreases. if > 0.4, doesnt work
        % this whole part must be written again, I think there is a
        % discontinuity in the profile. Maybe spline to smooth 
        beta0 = -( alpha_wanted.*r0^2)/(2*input.q(rSetAlpha)^2*input.R0*rSetAlpha);
        denomA = (1.1 - r0)^3 - (3/2)*(0.5-r0)*(1.1-r0)^2;
        aa = -2*beta0*(1.1-r0)/r0 - 2*beta0 + beta0*(1.1-r0)^2/(r0*(0.5-r0));
        aa = aa/denomA;
        bb = (-3*aa*(0.5 - r0)^2 - 2*beta0/r0)/(2*(0.5-r0));
        cc = beta0/(r0*bb);
        dd = 2*beta0 - bb*cc^2;
       
        betafter = @(x) aa.*(x-r0).^3 + bb.*(x-r0 + cc).^2+ dd;
        input.beta = @(x) (beta0*(1+(x./r0).^2)).*(x<=r0) + betafter(x).*(x>r0);
        input.betap = @(x) 2.*x.*beta0./r0^2.*(x<=r0) - 4*beta0.*(x-r0)./((1.05-r0).^2).*(x>r0);
        input.betapp = @(x) (2*beta0./r0^2).*(x<=r0)  - 4*beta0./((1.05-r0).^2).*(x>r0);
        
        input.alpha = @(x) -(input.q(x).^2).*input.R0.*input.betap(x);
        input.alphap = @(x) -2.*input.qp(x).*q(x).*input.R0.*input.betap(x) - input.q(x).^2.*input.R0.*input.betapp(x);
    end
end


% --- Handle optional parameters

if ~isempty(opts) && isfield(opts, 'additionalPlot') && toPlot
    input.additionalPlot = opts.additionalPlot;
else
    input.additionalPlot = 0; % to obtain plots of the q profile, pressure,
    %                           grid config.
end

%For m=1 and big alpha case, sometimes the main mode can couple a bit
% wierdly to the sideband. To avoid that, one can set
% opts.blockCouplingAfterRSp = 1 to decouple them right after the upper
% rational surface. Doesn't affect the growth rate.
if ~isempty(opts) && isfield(opts, 'blockCouplingAfterRSp')
    input.blockCouplingAfterRSp = opts.blockCouplingAfterRSp;
else
    input.blockCouplingAfterRSp = 0;
end

% Shaping of flux surfaces
if ~isempty(opts) && isfield(opts, 'addShaping')
    input.addShaping = opts.addShaping;
    if input.addShaping
        input.elong = profiles.elong;
        input.triangByR = profiles.triangByR;
    end
else
    input.addShaping = 0;
    input.elong = 1;
    input.triangByR = 0;
    warning('No shaping, if you want some, include opts.addShaping =1 (except for IM solvers)');
end

% If you want to retrieve modes associated with lower growth rates, hence 
% having more radial oscillations, set opts.higherOscModes = 1.
if ~isempty(opts) && isfield(opts,'higherOscModes')
    input.higherOscModes = opts.higherOscModes;
else
    input.higherOscModes = 0;
end

% Maximum number of iterations of Eigs solver before abandonning.
if ~isempty(opts) && isfield(opts, 'N_tries')
        input.N_tries = opts.N_tries;
else
        input.N_tries = 10;
end

% Choose boundary conditions for the upper sideband : Dirichlet (default),
% or Neumann if we stop the simulation just before the upper sideband
% rational surface. (This can be donne by setting opts.upperBound = 0.8,
% for instance, if you want to stop the simulation at r/a = 0.8)
input.NeumannBC = 0;
if ~isempty(opts) && isfield(opts, 'upperBound')
    input.upperBound = opts.upperBound;
    if isempty(input.RSp) || (input.upperBound <= RSp - 1e-4)  
        input.NeumannBC = 1;
        warning(['As the simulation upper boundary is within upper rational', ...
        'surface, Neumann BCs are assumed.']);
    end
else
    input.upperBound = 1;
end

input.reverse_qprofile = 0;
if ~isempty(opts) &&  isfield(opts, 'reverse_qprofile')
      if opts.reverse_qprofile == 1
         input.reverse_qprofile = 1;
      end
end

% Include Mercier (de)stabilisation on the sidebands
input.sidebandMercier = 0;
if ~isempty(opts) &&  isfield(opts, 'sidebandMercier')
    input.sidebandMercier = opts.sidebandMercier;
end


% --- Only needed for resistive simulations:

if resistiveSimu 

    % Lundquist number
    if (~isempty(profiles)) && isfield(profiles, 'SL') 
        input.SL = profiles.SL;
    elseif (~isempty(profiles)) && isfield(profiles, 'resistivity') 
        input.SL = mu0*RS^2*alfven(RS)/(profiles.resistivity);
    else
        input.SL = 10^7; 
    end
    
    input.resOnlyOnMainRatSurf = 0;
    input.resOnlyOnSidRatSurf = 0;

    if ~isempty(opts) && isfield(opts, 'resOnlyOnMainRatSurf')
        if opts.resOnlyOnMainRatSurf == 1
            input.resOnlyOnMainRatSurf = 1;
        end
    end
    if ~isempty(opts) &&  isfield(opts, 'resOnlyOnSidRatSurf')
          if opts.resOnlyOnSidRatSurf == 1
             input.resOnlyOnSidRatSurf = 1;
          end
    end

    % Model were we neglect perpendicular perturbation to potential vector
    % (roughly equals neglecting parallel magnetic field fluctuations)
    input.neglectDeltaBParallel = 0;
    if ~isempty(opts) &&  isfield(opts, 'neglectDeltaBParallel')
          if opts.neglectDeltaBParallel == 1
             input.neglectDeltaBParallel = 1;
          end
    end

    if isempty(input.RSp)
        input.resOnlyOnMainRatSurf = 1;
        warning(['No upper sideband rational surface in the plasma : resistivity only ...' ...
            'on the main rat. surf.'])
    end

    % Model where we do not include compression effects
    input.noCompr = 0;
    if ~isempty(opts) &&  isfield(opts, 'noCompr')
        if opts.noCompr == 1
            input.noCompr = 1;
        end
    end
end


end





