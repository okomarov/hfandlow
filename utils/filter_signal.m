function signal_filtered = filter_signal(mdl, signal)
[nrow,nser,nsig] = size(signal);
signal_filtered = NaN(nrow,nser,nsig);

for ii = 1:nsig
    ivalid = ~isnan(signal(:,:,ii));
    for c = 1:nser
        idx = ivalid(:,c);
        if nnz(idx) > 10
            y = signal(idx,c,ii);
            fit = mdl.estimate(y,'Display','off');
            signal_filtered(idx,c,ii) = y - fit.infer(y);
        end
    end
end
