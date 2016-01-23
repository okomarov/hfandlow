function out = iprice2dret(tb)
% IPRICE2DRET Convert intraday prices to daily returns
%
%   IPRICE2DRET(TB) TB is a table with serial 'Datetime' and 'Price'
%
%
%   OUT = ... 
%         OUT is a table with yyyymmdd 'Date' field and 'Ret'

if ~issorted(tb.Datetime)
    error('iprice2dret:unsorted','DATETIME should be sorted in ascending order.')
end

% Get rid of beginning of day NaNs
tb = tb(~isnan(tb.Price),:);

% Subscripts per day
[undates,~,subs] = unique(serial2yyyymmdd(tb.Datetime));

% Calculate returns
ret = accumarray(subs, tb.Price,[],@(x) x(end)/x(1)-1);

% Output
out = table(undates, ret, 'VariableNames',{'Date','Ret'});

end