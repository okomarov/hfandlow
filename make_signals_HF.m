function signals = make_signals_HF(permno,date,betacomp)

[unDt,~,midx] = unique(date/100);
nmonths       = numel(unDt);
nseries       = numel(permno);
signals       = NaN(nmonths, nseries,4);

% Monthly
for ii = 1:nmonths
    idx             = midx == ii;
    tmp             = squeeze(nansum(betacomp(idx,:,:)));
    betas           = tmp(:,1)./tmp(:,2);
    signals(ii,:,1) = 0.5*betas + 0.5*nanmean(betas);
end

% quarterly
for ii = 4:nmonths
    idx             = ismember(midx, ii-4+1:ii);
    tmp             = squeeze(nansum(betacomp(idx,:,:)));
    betas           = tmp(:,1)./tmp(:,2);
    signals(ii,:,2) = 0.5*betas + 0.5*nanmean(betas);
end

% semi-annually
for ii = 6:nmonths
    idx             = ismember(midx, ii-6+1:ii);
    tmp             = squeeze(nansum(betacomp(idx,:,:)));
    betas           = tmp(:,1)./tmp(:,2);
    signals(ii,:,3) = 0.5*betas + 0.5*nanmean(betas);
end

% Yearly
for ii = 12:nmonths
    idx             = ismember(midx, ii-12+1:ii);
    tmp             = squeeze(nansum(betacomp(idx,:,:)));
    betas           = tmp(:,1)./tmp(:,2);
    signals(ii,:,4) = 0.5*betas + 0.5*nanmean(betas);
end
end