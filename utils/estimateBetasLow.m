function betas = estimateBetasLow(lookback, rets, mkt)

% Unstack returns
warning off MATLAB:table:ModifiedVarnames
rets = rets(~isnan(rets.RetCC),{'Date','Permno','RetCC'});
rets = unstack(rets,'RetCC','Permno');
rets = sortrows(rets,'Date');
warning on MATLAB:table:ModifiedVarnames

% Intersect dates
[Date,ia,ib] = intersect(rets.Date, mkt.Date);
rets         = rets(ia,:);
mkt          = mkt(ib,:);

% Extract data
vnames = getVariableNames(rets);
rets   = rets{:,2:end};
rets2  = rets.^2;
mkt    = mkt{:,2:end};
mkt2   = mkt.^2;

[nobs,nseries] = size(rets);
betas          = NaN(nobs,nseries);
for row = lookback:nobs
    pos  = row-lookback+1:row;
    iobs = sum(~isnan(rets(pos,:))) > lookback * 0.5;
    
    % Covariance
    Exy = nanmean(bsxfun(@times, rets(pos,iobs), mkt(pos)));
    Ex  = nanmean(rets2(pos,iobs));
    Ey  = nanmean(mkt(pos));
    Cov = Exy - Ex*Ey;
    
    % Variance
    Ey2 = nanmean(mkt2(pos));
    Var = (Ey2 - Ey^2);
    
    % betas
    betas(row,iobs) = Cov./Var;
end
betas = nanfillts(betas);
betas = [table(Date), array2table(betas,'VariableNames',vnames(2:end))];

end