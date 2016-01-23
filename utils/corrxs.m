function corrmat = corrxs(panel, names)
% CORRXS Average of cross-sectional correlations
[nobs,~,nlay]  = size(panel);
if nargin < 2
    names = matlab.internal.table.dfltVarNames(1:nlay);
end
[pears, spear] = deal(zeros(nobs,nlay,nlay)); 
inan           = any(isnan(panel),3);
% inan(1:11,:)   = any(isnan(signals(1:11,:,1:end-1)),3); % First year mom12-2 is always NaN, skip
for ii = 1:nobs
    slice = squeeze(panel(ii,~inan(ii,:),:));
    if isempty(slice)
        continue
    end
    pears(ii,:,:) = corr(slice,'type','Pearson');
    spear(ii,:,:) = corr(slice,'type','Spearman');
end
pears   = squeeze(nanmean(pears));
spear   = squeeze(nanmean(spear));
corrmat = tril(pears,-1) + triu(spear,+1) + diag(NaN(nlay,1));
corrmat = array2table(corrmat, 'VariableNames', names, 'RowNames',names);
end