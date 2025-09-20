function[coefResMainChi, coefResSidChi] = MakeChiCoefNoSidRes(m, xall, RS, upperBound)
% If m=1, the ideal displacement stays non zero until r=0, and this creates
% problems, especially with the compressibility variable. So in this situation,
% we set to zero the resistive coefficient until just before the rational
% surface.

    if m==1
        lambda = 0.5.*tanh(-200.*(xall - (RS(1) - 0.05))) + 0.5.*tanh(200.*(xall - (upperBound - 0.01)));
    else
        lambda = 0.5.*tanh(-200.*xall) + 0.5.*tanh(200.*(xall - (upperBound - 0.01)));
    end
        coefResMainChi = exp(-10*lambda)/exp(10);
        coefResSidChi = 0.*xall;  
end
