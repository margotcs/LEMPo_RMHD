function xGrid = gridBunching(RS, RSp, RSm, N, upperBound, opts)
% Grid bunching around rational surfaces (main + sidebands) for FM models.

% INPUTS:
%   RS          - Main rational surface(s) (vector)
%   RSp         - Upper sideband rational surface(s) (scalar or vector)
%   RSm         - Lower sideband rational surface(s) (scalar or vector)
%   N           - Number of grid points
%   upperBound  - Radial end of the simulation (default is r/a = 1)
%   opts        - Option struct, to access linGrid option
%
% OUTPUT:
%   xGrid - grid points (column vector)


linGrid = 0;
if ~isempty(opts) && isfield(opts, 'linGrid') % if you just want a linear 
    % grid without bunching
    linGrid = opts.linGrid;
end


if ~linGrid

    bunchPoints = [double(RS); double(RSm); double(RSp)];
    bunchPoints = sort(bunchPoints);
    isBunchPointRS = bunchPoints; % this tab is going to be built just below
    NB = size(bunchPoints); NB = max(NB(:));
    
    % building isBunchPointRS tab, in which the entry will be 1 if the 
    % associated point is a main harmonic rational surface, ans zero 
    % otherwise. We need it because we bunch more around the main rational 
    % surface(s) than around the other locations.
    for j=1:NB % NB is the number of locations around which we bunch
        if min(abs(double(RS) - isBunchPointRS(j)))== 0
            isBunchPointRS(j) = 1; % means that this point is a main rational surface
        else
            isBunchPointRS(j) = 0;
        end
    end
    
    % We will bunch around NB + 1 points (all rational surfaces, plus a 
    % little bunch around r=0, so that eigenfunctions behave well there).
    % Thus, to get everything smooth we change N (number of points in the 
    % radial grid) a bit so that it is a multiple of NB + 1.
    NDividedByNB = round(N/(NB+1));   
    if N/(NB+1) ~= NDividedByNB
        N = NDividedByNB*(NB+1);
    end  
    %
    
    bTot = 1/0.00001;  % Variable to 'share' for the bunching. The lower it 
    % is, the most localized the bunching will be.
    NBMainRatSurf = nnz(isBunchPointRS);
    bunchZonesWidth = zeros(NB,1);
    cut = zeros(NB+1,1);
    cut(1) = bunchPoints(1)/2;
    
    % Divide the whole radial extent between the different zones where to
    % bunch.
    for i=1:NB
        if i ~= NB
            cut(i+1) = (bunchPoints(i+1) + bunchPoints(i))/2;
            bunchZonesWidth(i) = cut(i+1) - cut(i);
        else
            cut(i+1) = upperBound;
            bunchZonesWidth(i) = cut(i+1) - cut(i);
        end
    end
    % give to each points a factor of how much the grid will be bunched
    % aroung them, given the distance between each of them :
    bunchFactor = bunchZonesWidth*bTot;
    
    stock = 0; % stock variable is used to add more to the bunchFactor 
    % associated with the rational surface(s) of the MAIN harmonic. Like
    % this the bunching is more pronunced there.
    for i=1:NB
        if ~isBunchPointRS(i)
            bunchFactor(i) = bunchFactor(i)/2; 
            stock = stock + bunchFactor(i)/2; 
        end
    end
    
    for i =1:NB
        if isBunchPointRS(i)
            bunchFactor(i) = bunchFactor(i) + stock/NBMainRatSurf;
        end
    end
    
    xVec = zeros(NB, N/(NB+1));
    
    %first, bunching around zero:
    ppp = sym('ppp');
    b = 0.001; % !Do not change this value!
    l1 = vpasolve(ppp^3+b*ppp==0,ppp);
    l2 = vpasolve(ppp^3+b*ppp==cut(1),ppp);
    xTemp = linspace(l1(1),l2(1), N/(NB+1)); 
    xGrid = double(xTemp.^3 + b.*xTemp);
    step = xGrid(end) - xGrid(end-1);
    
    % Below : now we have everything we need, so this will be the actual
    % redefinition of the grid.
    for i=1:NB
        b = 1/bunchFactor(i);
            
        l1 = vpasolve(ppp^3+b*ppp + bunchPoints(i) ==cut(i) + step,ppp);
        l2 = vpasolve(ppp^3+b*ppp + bunchPoints(i) ==cut(i+1),ppp);
        xTemp = linspace(l1(1),l2(1), N/(NB+1)); 
        xVec(i,:) = double(xTemp.^3 + b.*xTemp + bunchPoints(i));
        step = xVec(i,end) - xVec(i,end-1);
    end
    
    xVec = reshape(xVec, 1, []);
    xVec = sort(xVec);
    xGrid = [xGrid xVec]';
else
    xGrid = linspace(0,upperBound,N)'; % simply a linear grid without bunching
end

end
