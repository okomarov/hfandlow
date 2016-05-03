function [iA, iS] = sample_refresh(A,S,minstep)
% SAMPLE_REFRESH Sample a pair of dates with the refresh-time scheme
%
%   SAMPLE_REFRESH(A,S)
%       A and S are vectors of serial dates
%
%   [iA, iS] = ...
%       The outputs are logical indices that select the sampled points
%
% Example:
%
%   A       = getTaqData('symbol','AAPL',20070906,20070906);
%   S       = getTaqData('symbol','SPY',20070906,20070906);
%   [ia,is] = sample_refresh(A.Datetime, S.Datetime, minstep);
%
%   References:
%   [1] Barndorff-Nielsen, O. E., Hansen, P. R., Lunde, A. & Shephard, N.
%       "Multivariate realised kernels: Consistent positive semi-definite
%        estimators of the covariation of equity prices with noise and
%        non-synchronous trading", Journal of Econometrics 162,
%        149–169 (2011)
%
% See also: DATENUM, FIXEDSAMPLING
if nargin < 3 || isempty(minstep)
    minstep = 0;
end

nA = numel(A);
nS = numel(S);
iA = false(nA,1);
iS = false(nS,1);
cA = 1;
cS = 1;
dateA = -1;
dateS = -1;

while cA <= nA && cS <= nS
    newDateA = A(cA);
    newDateS = S(cS);
    if newDateA - dateA < minstep && newDateS - dateS < minstep
        cA = cA + 1;
        cS = cS + 1;
        continue
    end

    dateA = newDateA;
    dateS = newDateS;

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
