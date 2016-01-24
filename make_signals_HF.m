function signals = make_signals_HF(permno,date,beta)

[unDt,~,midx] = unique(date/100);
nmonths       = numel(unDt);
nseries       = numel(permno);
signals       = NaN(nmonths, nseries,2);

% Monthly
for ii = 1:nmonths
    imonth          = midx == ii;
    tmp             = squeeze(nansum(beta(imonth,:,:)));
    betas           = tmp(:,1)./tmp(:,2);
    signals(ii,:,1) = 0.5*betas + 0.5*nanmean(betas);
end

% Yearly
for ii = 12:nmonths
    iyear           = ismember(midx, ii-12+1:ii);
    tmp             = squeeze(nansum(beta(iyear,:,:)));
    betas           = tmp(:,1)./tmp(:,2);
    signals(ii,:,2) = 0.5*betas + 0.5*nanmean(betas);
end
end