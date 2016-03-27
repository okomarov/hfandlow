function [coeff,se,tratio,pval] = regressHighOnLow(hfret,lfret)
ikeep = ~isnan(lfret);
nobs  = nnz(ikeep);
X     = [ones(nobs,1), lfret(ikeep)];
y     = hfret(ikeep);
opts  = {'bandwidth',floor(4*(nobs/100)^(2/9))+1,'intercept',false,'type','HAC','weights','BT','display','off'};

[~,se,coeff] = hac(X,y,opts{:});
tratio       = coeff./se;
pval         = 2 * normcdf(-abs(tratio));
end