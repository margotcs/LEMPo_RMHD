function makeFirstPlots(inputs, xall,N, resistiveSimu)
% Plots q profile, beta profile, alpha profile, and grid configuration.

% INPUTS:
%   inputs        - struct from getProfiles (contains q, beta, alpha, RS, etc.)
%   xall          - grid points (vector)
%   N             - number of grid points
%   resistiveSimu - boolean, to know if we included resistivity or not.



set(gca,'TickLabelInterpreter','latex')
set(0, 'DefaultAxesFontSize', 30)

% --- Plot q profile
figure
plot(xall, inputs.q(xall), 'LineWidth',3) ;
title(['q profile : m = ',num2str(inputs.m),' ,n= ',num2str(inputs.n)],'interpreter',...
        'tex', 'FontWeight', 'normal');
xline(inputs.RS, '--b', 'LineWidth', 2);
if ~isempty(inputs.RSp)
    xline(inputs.RSp, '--r', 'LineWidth', 2);
end
xlabel('r/a');

% --- Plot grid configuration
figure
plot( linspace(0,1,N), xall, 'LineWidth',3);
title('Grid configuration')

% --- Plot beta profile
figure
plot( xall, inputs.beta(xall),'LineWidth',3);
title('\beta profile')

% --- Plot alpha profile

figure
plot( xall, inputs.alpha(xall),'LineWidth',3);
hold on
plot(xall, inputs.alpha(xall)./xall,'--','LineWidth',3);
legend('\alpha', '\alpha / r')
title('\alpha profile')

if resistiveSimu
    cl = EPFLcolors();  % Get the struct with all colors
    colortab = {'groseille', 'canard', 'rose', 'zinzolin', 'chartreuse', 'montRose', 'ardoise', 'vertDeau'};
    colormat = cell2mat(cellfun(@(c) EPFLcolors.(c), colortab, 'UniformOutput', false));


    figure
    set(gca,'ColorOrder',EPFLcolors(colortab),'nextplot','replacechildren')
    set(gca,'FontName','Cambria Math' ) 

    plot(xall, coefResMainChi, xall, coefResSidChi, 'LineWidth',3);
    for i=1:max(size(RS))
       xline(double(RS(i)), '--', 'LineWidth', 2, 'Color','#B51F1F');
    end
    if ~isempty(RSp)
        xline(RSp, '--', 'LineWidth', 2,'Color', '#007480');
    end
    if ~isempty(RSm)
        xline(RSm, '--', 'LineWidth', 2);
    end
    legend('Coef. for  \chi^{(m)} in equations', 'Coef. for  \chi^{(m+1)} in equations',...
            'r_m', 'r_{m+1}', 'r_{m-1}');
    xlabel('r/a')
end


end