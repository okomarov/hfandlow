function [priceA, priceS] = sample_preaverage(dtA,dtS,priceA,priceS,minstep)
% SAMPLE_PREAVERAGE Pre-average price series at refresh-time sampling points
%
%   SAMPLE_PREAVERAGE(DTA,DTS,PRICEA,PRICES)
%       DTA and DTS are vectors of serial dates and PRICEA and PRICES are
%       vectors of prices
%
%   [PRICEA, PRICES] = ...
%       The outputs are vectors of pre-averaged prices at refresh-time
%       sampling points
%
% Example:
%
%   A = getTaqData('symbol','AAPL',20070906,20070906);
%   S = getTaqData('symbol','SPY',20070906,20070906);
%   [priceA,priceS] = sample_preaverage(A.Datetime,S.Datetime,A.Price,S.Price,minstep);
%
%
%   References:
%   [1] K. Christensen, S. Kinnebrock, & M. Podolskij, "Pre-averaging 
%        estimators of the ex-post covariance matrix in noisy diffusion 
%        models with non-synchronous data", Journal of Econometrics 159, 
%        116-133 (2010)
%   [2] Barndorff-Nielsen, O. E., Hansen, P. R., Lunde, A. & Shephard, N.
%       "Multivariate realised kernels: Consistent positive semi-definite
%        estimators of the covariation of equity prices with noise and
%        non-synchronous trading", Journal of Econometrics 162,
%        149–169 (2011)
%
% See also: SAMPLE_REFRESH, FIXEDSAMPLING

% Refresh points
[ia,is] = sample_refresh(dtA, dtS, minstep);
priceA  = priceA(ia);
priceS  = priceS(is);

% Last point to have k prices ahead
n   = numel(priceA);
k   = floor(sqrt(n));
pos = (1:n-k+1)';

% Average
pos    = bsxfun(@plus,pos,0:k-1);
priceA = sum(priceA(pos),2)/k;
priceS = sum(priceS(pos),2)/k;
end
