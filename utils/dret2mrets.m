function [mrets, mdates] = dret2mrets(dates,ret)
if ~isdatetime(dates)
    dates = yyyymmdd2datetime(dates);
end
sz     = size(ret);
idx    = [true; logical(diff(month(dates)))];
rsub   = repmat(cumsum(idx), 1, sz(2));
csub   = repmat(1:sz(2),sz(1),1);
mrets  = accumarray([rsub(:),csub(:)], ret(:), [], @(r) prod((1+r))-1);
mdates = dates(idx);
end