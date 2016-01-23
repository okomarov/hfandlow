function rs = runsum(lookback, data)
rs        = NaN(size(data));
nonan     = ~isnan(data);
rs(nonan) = filter(ones(lookback,1), 1, data(nonan), NaN(lookback-1,1));
rs        = {rs};
end