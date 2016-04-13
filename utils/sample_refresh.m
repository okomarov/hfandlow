function [iA, iS] = sample_refresh(A,S)
% SAMPLE_REFRESH Sample a pair of dates with the refresh-time scheme
%
%   SAMPLE_REFRESH(A,S)
%       A and S are vectors of serial dates
%
%   [iA, iS] = ...
%       The outputs are logical indices that select the sampled points
%
%
%   References:
%   [1] Barndorff-Nielsen, O. E., Hansen, P. R., Lunde, A. & Shephard, N.
%       "Multivariate realised kernels: Consistent positive semi-definite
%        estimators of the covariation of equity prices with noise and
%        non-synchronous trading", Journal of Econometrics 162,
%        149–169 (2011)
%
% See also: DATENUM, FIXEDSAMPLING

nA = numel(A);
nS = numel(S);
iA = false(nA,1);
iS = false(nS,1);
cA = 1;
cS = 1;

while cA < nA && cS < nS
    dateA = A(cA);
    dateS = S(cS);

    if dateA > dateS
        iA(cA) = true;
        try
            while S(cS) <= dateA
                cS = cS + 1;
            end
        catch
        end
        iS(cS-1) = true;

        cA = cA+1;

    else
        iS(cS) = true;
        try
            while A(cA) <= dateS
                cA = cA + 1;
            end
        catch
        end
        iA(cA-1) = true;

        cS = cS+1;
    end
end
end
