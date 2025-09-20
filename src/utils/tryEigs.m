function [V, ev] = tryEigs(A, B, ev_guess, maxNbTry, xall, N, rStar, reverse_qprofile, ...
    higherOscModes, additionalPlot)
% Attempts (maxNbTry times) to solve the generalised eigenvalue proble 
% B*V*ev = A*V 

% INPUTS:
%   A, B                - Sparse Matrices defining the problem
%   ev_guess            - Initial guess for the eigenvalue (should be much
%                         above the expected value)
%   maxNbTry            - Maximum number of running eigs before stopping
%   xall                - Grid
%   N                   - Number of points in the grid
%   rStar               - Characteristic length (1st rational surface),
%                         useful for some checks on the eigenfunction
%   reverse_qprofile    - Boolean: 1 if the q profile is reversed (useful
%                         for checks)
%   higherOscModes      - Boolean: 1 if you want to retrieve higher
%                         oscillation (less unstable) modes. If 0 we look 
%                         for the most unstable mode
%   additionalPlot      - Boolean: 1 if you want additional plots of what
%                         happened during the attempts to use eigs.
%
% OUTPUTS:
%   V    - Eigenvectors (radial plasma displacements of main harmonic and
%          sidebands, resistive and compressibility corrections, and dummy
%          variables)
%   ev   - Eigenvalue (will be growth rate for resistive cases, and growt
%          rates square for ideal cases)

notconverge =1;
count = 0; % counting the number of iterations before stopping
result.warning = 0;

while notconverge
    count = count + 1;
    if count > maxNbTry
        warning('Too many tries, guess %e is not convenient. See previous displays to modify it. Output will be ev=0.', ev_guess);
        V = [];
        ev = 0;
        break;
    end

    try
        [V,D] = eigs(A ,B ,1,ev_guess);       
        ev = diag(D)
        [ev,I]=sort(real(ev),'descend');
        ev = ev(1);
        V = V(:,I(1));
    catch
        D = -0.01 ; 
        ev = -0.01;
    end


   if (ev == -0.01) || isnan(ev) || ~isreal(D)
        warning('Guess %e was too high, new guess is :', ev_guess);
        ev_guess = 0.7*ev_guess
        result.warning = 1;
   else
       Xi = V(1:N);
        [maxXi, idx] = max(abs(Xi));
        if Xi(idx) < 0
            Xi = - Xi;
        end
        Xi = Xi./maxXi;

        % check that we dont obtain the non physical solution where \xi goes
        % first negative : ( will cross zero once and only once before r_s )
        [~, idxRS] = min(abs(xall - rStar));
        [reszero, ~] = min(abs(Xi(round(idxRS/2):idxRS)));
        
        if min(Xi(1:idxRS)) < -0.01
            goesbelowzero = 1;
        else
            goesbelowzero = 0;
        end
         % check that it doesnt stay equal to zero all the time 
         % we don't perform this check if we have a reverse q profile.
         % indeed often the first rat. surf will be too close to the
         % axis for the check to be meaningful. 
         if (reverse_qprofile) == 0 && ( (max(abs(diff(Xi(round(idxRS/2):idxRS))./diff(xall(round(idxRS/2):idxRS)))) <= 0.05 ) || (max(abs(Xi(round(idxRS/2):idxRS))) <= 0.1))
               isFlat = 1;
         else
               isFlat =0;
         end
        

         if ~higherOscModes && (goesbelowzero && ~isFlat) 
            % Here there are two different cases, first one is the one 
            % having a lot of oscillations, second one going below zero 
            % before having a kind of normal behaviour. Both are non 
            % physical so we want to discard them.
            [pks, locs] = findpeaks(Xi,xall,'MinPeakHeight',0.05);
            sz = size(locs); sz = max(sz);
            if sz > 5
                result.warning =1;
                warning('Guess %e was too high. New guess is :', ev_guess);
                if additionalPlot
                    figure
                    plot(xall, Xi)
                    title(['\xi, the guess ', num2str(ev_guess),' was too high.'])
                end   
                ev_guess = 0.7*ev_guess
             else
                result.warning =1;
                warning('Guess %e was too low. New guess is :', ev_guess);
                if additionalPlot
                    figure
                    plot(xall, Xi)
                    title(['\xi, the guess ', num2str(ev_guess),' was too low.'])
                end   
                ev_guess = 1.5*ev_guess
              end
         elseif (reszero < 0.01) && isFlat
             warning('Guess %e was too high, new guess is :', ev_guess);
             result.warning = 1;
             if additionalPlot
                figure
                plot(xall, Xi)
                title(['\xi, the guess ', num2str(ev_guess),' was too high.'])
             end 
            ev_guess = 0.7*ev_guess
         else 
            notconverge = 0;
            if additionalPlot
                figure
                plot(xall, Xi)
                title('\xi, guess is good')
            end 
        end
    end
end   

end