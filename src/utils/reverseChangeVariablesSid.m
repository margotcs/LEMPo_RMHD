function XiFinal = reverseChangeVariablesSid(xall, Xi, Chi, XiBar, N, ...
    sidebandKind,ev, m, dqqs, deltap, rStar, RSp, NeumannBC, resistiveSimu,...
    additionalPlot)

% We computed everything for Xip bar, so variable Xip is actually
% the bar variable for now (Respectively Xim bar, and Xim for the lower
% sideband). So this function returns to the original sidebands expression. 
% (See below eq. 10 in paper "Fundamental properties of ideal and resistive 
% infernal modes in tokamaks"

% INPUTS:
%   xall                - Grid
%   Xi                  - Main harmonic radial plasma displacement
%   Chi                 - Resistive variable (if ideal simu, set to [])
%   XiBar               - Sideband output from eigs: so it's Xim bar, or
%                         Xip bar
%   N                   - Number of points in the grid
%   sidebandKind        - String: either "upper" (m+1 sideband), or "lower"
%                         (m-1 sideband)
%   ev                  - \gamma / \omega_A associated with the mode
%   m                   - Poloidal mode number
%   dqqs                - Function handle for \Delta q/qs
%   deltap              - Function handle for derivatives of Shafranov
%                         shift
%   rStar               - Characteristic length (scalar, 1st main rational
%                         surface)
%   RSp                 - Upper sideband rational surface (can be [])
%   NeumannBC           - Boolean: 1 if we stop before upper sideband
%                         rational surface and apply Neumann BCs there
%   resistiveSimu       - Boolean: 1 if we are running a resistive model
%   additionalPlot      - Boolean: 1 if you want plot of the change of
%                         variable process for the sidebands
%
% OUTPUTS:
%   XiFinal             - Actual sideband radial displacement


if resistiveSimu
    xiR_fun = dqqs(xall).*Xi - coefResMainChi.*Chi; % resistive case
else
    xiR_fun = dqqs(xall).*Xi; % ideal case
end

 xiRp = diff(xiR_fun)./diff(xall);
 xiRp = [0; xiRp];
 
 switch sidebandKind
     case "upper"
         intXiR = (1+m).*cumtrapz(xall, xall.^(2+m).*deltap(xall).*xiRp);
         Xi_int = XiBar + intXiR./xall.^(2+m); %integral from the left
     case "lower"
         intXiR = (1-m).*cumtrapz(flip(xall), flip(xall.^(2-m).*deltap(xall).*xiRp));
         intXiR = flip(intXiR)./(xall.^(2-m));
         Xi_int = XiBar - intXiR; %integral from the right
         XiFinal = Xi_int;
 end


if strcmp(sidebandKind, "upper")
 % for the upper sideband sometimes problem because of the integral in the
 % change of variable definition, so need to do additional steps.
    if (m==1) || (m==2) || (ev > 0.1)
        if NeumannBC || ~exist('RSp','var') || isempty(RSp)
            [~, indRS] = min(abs(xall - rStar ));
            [~, indmatchPlus] = min(abs(Xi_int((indRS + round(N/20)):end) - XiBar((indRS + round(N/20)):end)));  
            indmatchPlus = indRS + round(N/20) + indmatchPlus;
        else
            [~, indRSp] = min(abs(xall - RSp));
            [~, indmatchPlus] = min(abs(Xi_int((indRSp + round(N/20)):end) - XiBar((indRSp + round(N/20)):end)));  
            indmatchPlus = indRSp + round(N/20) + indmatchPlus;
        end   
    if (indmatchPlus >= N-5) % in case there is no problem with the BC (lower growth rates )
        XiFinal = Xi_int;
    else
        XiFinal = [Xi_int(1:indmatchPlus); XiBar(indmatchPlus+1:N)];
    end
   else
     XiFinal = Xi_int;  % case where there is no need for additional manip.
   end
end   
     
if additionalPlot
     figure
     plot(xall, XiBar, '-', xall, XiFinal, '--', xall, Xi_int, '--','LineWidth', 3);
     legend('$\overline{\xi^{m \pm 1}}$', '$\xi^{m \pm 1}$', '\xi int')
end


end
