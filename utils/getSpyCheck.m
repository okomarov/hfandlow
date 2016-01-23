function getSpyCheck(freq)

% Take last price from spy sampled
spy = getSpy(freq);
spy = sortrows(spy,'Datetime');
pos = [find(diff(int64(serial2yyyymmdd(spy.Datetime))) ~= 0); size(spy,1)];
spy = spy(pos,:);

% Load dsfquery daily prices
dsfquery = loadresults('dsfquery');
dsfquery = dsfquery(dsfquery.Permno == 84398,:);
dsfquery = sortrows(dsfquery, 'Date');

% Intersect dates
[~,ia,ib] = intersect(serial2yyyymmdd(spy.Datetime), dsfquery.Date);
spy       = spy(ia,:);
dsfquery  = dsfquery(ib,:);

% Visual inspection
subplot(211)
plot(serial2datetime(spy.Datetime), spy.Price,...
     yyyymmdd2datetime(dsfquery.Date), abs(dsfquery.Prc))
legend({'Sampled','Dsfquery'})

p2r    = @(x) x(2:end)./x(1:end-1)-1;
perATE = abs(p2r(spy.Price) - p2r(abs(dsfquery.Prc)));
subplot(212);
boxplot(perATE)
str    = {'Absolute Tracking Error'
       sprintf('%-10s%5.4f','mean:',mean(perATE))
       sprintf('%-10s%5.4f','std:',std(perATE))};
title(str)

% Draw unsampled random day
date                       = serial2yyyymmdd(randsample(spy.Datetime,1));
unsampled                  = getSpy(inf,date,date);
unsampled                  = unsampled(~selecttrades(unsampled),:);
[Datetime,~,subs]          = unique(unsampled.Datetime);
Price                      = accumarray(subs, unsampled.Price,[],@fast_median);
% Sample with MFE
pt                         = 'D:\TAQ\HFbetas\utils\MFE'; addpath(pt) 
grid                       = (9.5/24:freq/(60*24):16/24)';
[MfePrice, ~, MfeDatetime] = realized_price_filter(double(Price), Datetime,'unit','fixed', grid+yyyymmdd2serial(date));

% Same sampled day
sampled = getSpy(freq,date,date);

% Compare
tb                           = [sampled, table(MfeDatetime,MfePrice)];
tb.MfePrice(isnan(tb.Price)) = NaN;

plot(serial2datetime(tb.Datetime),tb{:,{'Price','MfePrice'}})
end