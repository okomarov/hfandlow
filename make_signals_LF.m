function [signals,hpr,rf,mdate] = make_signals_LF(ret,date,factors)
[mdate,~,midx] = unique(date/100);
nmonths       = numel(mdate);
nseries       = size(ret,2);

signals = NaN(nmonths, nseries,2);
hpr     = NaN(nmonths, nseries);
rf      = NaN(nmonths, 1);

% Monthly
for ii = 1:nmonths
    imonth = midx == ii;
    nobs   = nnz(imonth);

    % Month of all non-nan returns
    r      = ret(imonth,:);
    nonans = all(~isnan(r));
    ngood  = nnz(nonans);
    r      = r(:,nonans);
    
    % Betas
    Y     = bsxfun(@minus,r, factors.RF(imonth)/100);
    X     = [ones(nobs,1), factors.MktMinusRF(imonth)/100];
    coeff = NaN(2,ngood);

    for jj = 1:ngood
        coeff(:,jj) = X\Y(:,jj);
    end
    betas = 0.5*coeff(2,:) + 0.5*mean(coeff(2,:));
    signals(ii,nonans,1) = betas;
end

% Yearly
for ii = 12:nmonths
    iyear = ismember(midx, ii-12+1:ii);
    nobs  = nnz(iyear);

    % Month of all non-nan returns
    r      = ret(iyear,:);
    inan   = isnan(r);
    nonans = sum(~inan) >= 200;
    ngood  = nnz(nonans);
    inan   = inan(:,nonans);
    r      = r(:,nonans);

    % Betas
    Y     = bsxfun(@minus,r, factors.RF(iyear)/100);
    X     = [ones(nobs,1), factors.MktMinusRF(iyear)/100];
    coeff = NaN(2,ngood);

    for jj = 1:ngood
        idx         = ~inan(:,jj);
        coeff(:,jj) = X(idx,:)\Y(idx,jj);
    end
    
    % Shrink towards xs mean
    betas = 0.5*coeff(2,:) + 0.5*mean(coeff(2,:)); 
    signals(ii,nonans,2) = betas;
end

% Holding period return
for ii = 1:nmonths
    imonth = midx == ii;

    r                 = ret(imonth,:);
    inan              = isnan(r);
    r(inan)           = 0;
    hpr(ii,:)         = prod(1+r)-1;
    hpr(ii,all(inan)) = NaN;

    rf(ii) = prod(1+factors.RF(imonth)/100)-1;
end
end