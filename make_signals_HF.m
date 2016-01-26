function signals = make_signals_HF(permno,date,betacomp,w)

[unDt,~,midx] = unique(date/100);
nmonths       = numel(unDt);
nseries       = numel(permno);
signals       = NaN(nmonths, nseries,4);

signals(:,:,1) = getBetaHF(nmonths,nseries,betacomp,midx,1,w(1));
signals(:,:,2) = getBetaHF(nmonths,nseries,betacomp,midx,4,w(2));
signals(:,:,3) = getBetaHF(nmonths,nseries,betacomp,midx,6,w(3));
signals(:,:,4) = getBetaHF(nmonths,nseries,betacomp,midx,12,w(4));
end

function betas = getBetaHF(nmonths,nseries,betacomp,midx,len,w)
betas = NaN(nmonths, nseries);
for ii = len:nmonths
    idx         = ismember(midx, ii-len+1:ii);
    tmp         = squeeze(nansum(betacomp(idx,:,:)));
    betas(ii,:) = tmp(:,1)./tmp(:,2);
end
if w ~= 1
    betas = bsxfun(@plus,w*betas, (1-w)*nanmean(betas,2));
end
end