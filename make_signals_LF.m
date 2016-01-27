function [signals,hpr,rf,mdate] = make_signals_LF(ret,date,factors,w,lag)
if nargin < 4
    isBayesW = true;
else
    isBayesW = false;
end
if nargin < 5 || isempty(lag)
    lag = 0;
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
    signals(ii,:,1) = getBetasLF(ret, factors, idx, lag, w(1), isBayesW);

    % HPR
    r                 = ret(idx,:);
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
    signals(ii,:,2) = getBetasLF(ret, factors, idx, lag, w(2), isBayesW, 50);
end
% Semi-annual
for ii = 6:nmonths
    idx             = ismember(midx, ii-6+1:ii);
    signals(ii,:,3) = getBetasLF(ret, factors, idx, lag, w(3), isBayesW, 100);
end
% Annual
for ii = 12:nmonths
    idx             = ismember(midx, ii-12+1:ii);
    signals(ii,:,4) = getBetasLF(ret, factors, idx, lag, w(4), isBayesW, 200);
end
end

function betas = getBetasLF(ret, factors, idx, lag, w, isBayesW, minobs)
r     = ret(idx,:);
ff    = factors(idx,:);
betas = NaN(1+lag,size(r,2));
n     = size(r,1);

% Dimson betas if lag > 1
for jj = 0:lag
    if nargin < 7
        minobs = n-jj;
    end
    [betas(jj+1,:),sigma_ts] = getBetasLF_(r(1+jj:n,:),ff(1:n-jj,:),minobs,isBayesW);
end
betas = sum(betas,1);

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
Y        = bsxfun(@minus,r, factors.RF/100);
X        = [ones(size(r,1),1), factors.MktMinusRF/100];
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