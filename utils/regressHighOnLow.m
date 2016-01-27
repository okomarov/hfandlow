function [coeff,se,tratio,pval] = regressHighOnLow(hcell,lcell,idx)
sz = size(hcell{1});
if nargin < 3
    idx = true(sz(1),1);
end
opts = {'intercept',false,'type','HAC','weights','BT','display','off'};
nsig = numel(hcell);
l    = ones(sz(1),1);

[coeff,se] = deal(NaN(nsig,2));
for ii = 1:nsig
    isel = ~isnan(lcell{ii}(:,end)) & idx(:);
    nobs = nnz(isel);
    X    = [l(isel), lcell{ii}(isel,end)];
    y    = hcell{ii}(isel,end);
    
    [~,se(ii,:),coeff(ii,:)] = hac(X,y,opts{:}, 'bandwidth',floor(4*(nobs/100)^(2/9))+1);
end
% tstat and pvalues
tratio = coeff./se;
pval   = 2 * normcdf(-abs(tratio));
end