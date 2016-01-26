function signals = make_signals_HF(permno,date,betacomp,w)

[unDt,~,midx] = unique(date/100);
nmonths       = numel(unDt);
nseries       = numel(permno);
signals       = NaN(nmonths, nseries,4);

signals(:,:,1) = getBetaHF(nmonths,nseries,betacomp,midx,1,w);
signals(:,:,2) = getBetaHF(nmonths,nseries,betacomp,midx,4,w);
signals(:,:,3) = getBetaHF(nmonths,nseries,betacomp,midx,6,w);
signals(:,:,4) = getBetaHF(nmonths,nseries,betacomp,midx,12,w);
end

function betas = getBetaHF(nmonths,nseries,betacomp,midx,len,w)
betas = NaN(nmonths, nseries);
for ii = len:nmonths
    idx         = ismember(midx, ii-len+1:ii);
    tmp         = squeeze(nansum(betacomp(idx,:,:)));
    betas(ii,:) = tmp(:,1)./tmp(:,2);
end
if w ~= 1
    betas = w*betas + (1-w)*nanmean(betas,2);
end
end