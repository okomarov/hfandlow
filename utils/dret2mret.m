function [mret, mdts] = dret2mret(ret, dts, refdts, requireall, endskip, startskip)
% [mret, mdts] = dret2mret(ret, dts, refdts)

sz = size(ret);
if numel(dts) ~= sz(1)
    error('d2mret:numObs','RET must have same number of rows as dates in DTS.')
end
if numel(sz) > 2
    colons = repmat({':'}, sz(3),1);
else
    colons = {':'};
end

% Intersect dates with reference dates at the monthly freq
if nargin >= 3 && ~isempty(refdts)
    idx = ismember(dts/100,refdts/100);
    dts = dts(idx);
    ret = ret(idx,colons{:});
end

% Monthly indexing with last day of month
[~,lastpos,msubs] = unique(dts/100,'last');
mdts              = dts(lastpos);
firstpos          = [1; lastpos(1:end-1)+1];

% Preallocate
N    = numel(lastpos);
mret = NaN([N,sz(2:end)]);

% Convert all NaNs to zero
inan = isnan(ret);
if ~requireall
    ret = nan2zero(ret);
end

% Monthly returns
for ii = 1:N
    % Select current month
    imonth      = ii == msubs;
    pos         = lastpos(ii)-endskip+1:lastpos(ii);
    imonth(pos) = 0;
    pos         = firstpos(ii):firstpos(ii)+startskip-1;
    imonth(pos) = 0;
    mret(ii,colons{:})  = prod(ret(imonth,colons{:}) + 1,1)-1;
    if ~requireall
        inanm          = all(inan(imonth,colons{:}),1);
        mret(ii,inanm) = NaN;
    end
end
end