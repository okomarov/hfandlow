function [signals,hpr,rf,mdate] = make_signals_LF(ret,date,factors)
[mdate,~,midx] = unique(date/100);
nmonths        = numel(mdate);
nseries        = size(ret,2);

signals = NaN(nmonths, nseries,4);
hpr     = NaN(nmonths, nseries);
rf      = NaN(nmonths, 1);

% Monthly
for ii = 1:nmonths
    idx             = midx == ii;
    signals(ii,:,1) = getBetasLF(ret(idx,:),factors(idx,:),nnz(idx));

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
    signals(ii,:,2) = getBetasLF(ret(idx,:),factors(idx,:),50);
end
% Semi-annual
for ii = 6:nmonths
    idx             = ismember(midx, ii-6+1:ii);
    signals(ii,:,3) = getBetasLF(ret(idx,:),factors(idx,:),100);
end
% Annual
for ii = 12:nmonths
    idx             = ismember(midx, ii-12+1:ii);
    signals(ii,:,4) = getBetasLF(ret(idx,:),factors(idx,:),200);
end
end

function betas = getBetasLF(r,factors,minobs)
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
coeff    = NaN(1,ngood);
sigma_ts = NaN(1,ngood);

for jj = 1:ngood
    idx          = ~inan(:,jj);
    n            = nnz(idx);
    x            = X(idx,:);
    y            = Y(idx,jj);
    b            = x\y;
    res          = y-x*b;
    Sxx          = sum((x(:,2)-sum(x(:,2))/n).^2);
    sigma_ts(jj) = sum(res.^2)/(n-2)/Sxx;
    coeff(jj)    = b(2);
end
sigma_xs = var(coeff);
w        = 1-sigma_ts./(sigma_ts+sigma_xs);

% Shrink towards xs mean
betas(nonans) = w.*coeff + (1-w)*mean(coeff);
end