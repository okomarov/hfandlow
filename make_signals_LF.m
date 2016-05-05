function [signals,hpr,rf,mdate] = make_signals_LF(ret,date,factors,w)
if nargin < 4
    isBayesW = true;
else
    isBayesW = false;
end

[mdate,~,midx] = unique(date/100);
nmonths        = numel(mdate);
nseries        = size(ret,2);

signals = NaN(nmonths, nseries,4);
hpr     = NaN(nmonths, nseries);
rf      = NaN(nmonths, 1);

% Monthly
for ii = 1:nmonths
    idx             = midx == ii;
    r               = ret(idx,:);
    signals(ii,:,1) = getBetasLF(r, factors(idx,:), w(1), isBayesW, size(r,1));

    % HPR
    inan              = isnan(r);
    r(inan)           = 0;
    hpr(ii,:)         = prod(1+r)-1;
    hpr(ii,all(inan)) = NaN;

    % RF
    rf(ii) = prod(1+factors.RF(idx)/100)-1;
end
% Quarterly
for ii = 4:nmonths
    idx             = ismember(midx, ii-4+1:ii);
    signals(ii,:,2) = getBetasLF(ret(idx,:), factors(idx,:), w(2), isBayesW, 50);
end
% Semi-annual
for ii = 6:nmonths
    idx             = ismember(midx, ii-6+1:ii);
    signals(ii,:,3) = getBetasLF(ret(idx,:), factors(idx,:), w(3), isBayesW, 100);
end
% Annual
for ii = 12:nmonths
    idx             = ismember(midx, ii-12+1:ii);
    signals(ii,:,4) = getBetasLF(ret(idx,:), factors(idx,:), w(4), isBayesW, 200);
end
end

function betas = getBetasLF(r, ff, w, isBayesW, minobs)

[betas,sigma_ts] = getBetasLF_(r, ff, minobs, isBayesW);

if isBayesW
    sigma_xs = nanvar(betas,2);
    w        = 1-sigma_ts./(sigma_ts+sigma_xs);
end

if w ~= 1
    betas = w.*betas + (1-w);
%     betas = w.*betas + (1-w)*nanmean(betas,2);
end
end

function [betas,sigma_ts] = getBetasLF_(r,factors,minobs,isBayesW)
betas = NaN(1,size(r,2));

% At least x non NaN returns
inan   = isnan(r);
nonans = sum(~inan) >= minobs;
ngood  = nnz(nonans);
inan   = inan(:,nonans);
r      = r(:,nonans);

% Betas
Y        = log(bsxfun(@minus,r, factors.RF/100)+1);
X        = [ones(size(r,1),1), log(factors.MktMinusRF/100+1)];
inanmkt  = isnan(factors.MktMinusRF);
coeff    = NaN(1,ngood);
sigma_ts = NaN(1,ngood);

for jj = 1:ngood
    idx = ~(inan(:,jj) | inanmkt);
    n   = nnz(idx);
    x   = X(idx,:);
    y   = Y(idx,jj);
    b   = x\y;
    if isBayesW
        res          = y-x*b;
        Sxx          = sum((x(:,2)-sum(x(:,2))/n).^2);
        sigma_ts(jj) = sum(res.^2)/(n-2)/Sxx;
    end
    coeff(jj) = b(2);
end
betas(nonans) = coeff;
end