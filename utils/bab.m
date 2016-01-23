function [ptfret,z, wl, wh, avgbeta] = bab(ret, betas, rf)
% BAB Betting-against-beta portfolio returns
%
%   BAB(RET, BETAS, RF) 
%       Same size RET and BETAS panels with the first row of RET
%       ahead in time, e.g. t+1, and the first row in BETAS as 
%       current time, i.e. t.
%       The risk-free, RF, should have same number of observations 
%       and same time-indexing as RET.
%       
%
%   [PTFRET, Z, WL, WH, AVGSIG] = ... 
%       PTFRET  - returns of the [low-beta, high-beta]
%       Z       - beta ranks by row
%       WL/WH   - weight matrices same size as BETAS
%       AVGBETA - returns the series of average [low-beta, high-beta]
%
%   Follows methodology for one-period, hence time subscript is suppressed.
%
%   With information up to time t we have and betas_t a 1xk vector:
%       z    = rank(betas_t)
%       w_h  =  max(z-zbar,0)/k
%       w_l  = -min(z-zbar,0)/k
%       k    = sum(|z-zbar|,2)/2
%   Then the return on the BAB portfolio, in the next period, is:
%       retBAB_t+1 =   (ret_t+1-rf)*w_l' / (betas_t*w_l')
%                    - (ret_t+1-rf)*w_h' / (betas_t*w_h')
%
% Reference:
%   2014 Frazzini, Pedersen - Betting against beta - JFE

% Bab
[nobs,nseries] = size(betas);

% rank
z = NaN(nobs,nseries);
for ii = 1:nobs
    z(ii,:) = tiedrank(betas(ii,:));
end

% Calculate zbar, z-zbar and normalizing k
zbar       = (max(z,[],2)+1)/2;
inan       = isnan(zbar);
zMinusZbar = bsxfun(@minus,z,zbar);
k          = nansum(abs(zMinusZbar)/2,2);

wh = bsxfun(@rdivide,  max(zMinusZbar,0), k);
wl = bsxfun(@rdivide, -min(zMinusZbar,0), k);

rh     = nansum(ret  .* wh,2);
bh     = nansum(betas.* wh,2);
rl     = nansum(ret  .* wl,2);
bl     = nansum(betas.* wl,2);
ptfret = [(rl-rf)./bl (rh-rf)./bh];

ptfret(inan,:) = NaN;

if nargout == 5
    avgbeta = [bl, bh];
    avgbeta(inan,:) = NaN;
end
end
