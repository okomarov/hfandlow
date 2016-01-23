%%

addpath 'D:\TAQ\HF\hfandlow\utils\MFE'

fprintf('Checking one random beta.\n')
beta   = getBetas(1,5,true,false,false,true,true);
record = beta(randsample(size(beta,1),1),:);

date   = record.Date;
permno = record.Permno;
% record = beta(beta.Permno == permno & beta.Date == date,:);

grid = (9.5/24:5/(60*24):16/24)';
ret  = loadresults('return_intraday_overnight');

%%%%%%%%%%%%%%%%%%
% Num - covariance
%%%%%%%%%%%%%%%%%%

% Unsampled series
unsampled         = getTaqData('permno',permno,date,date);
unsampled         = unsampled(~isInvalidTrade(unsampled),:);
[Datetime,~,subs] = unique(unsampled.Datetime);
Price             = accumarray(subs, unsampled.Price,[],@fast_median);

% Sample with MFE
Price = realized_price_filter(double(Price), Datetime,...
                              'unit','fixed', grid+yyyymmdd2serial(date));

% Inherit same NaNs
s           = getTaqData('permno',permno,date,date,[],'D:\TAQ\HF\data\TAQ\sampled\5min\nobad_vw');
inan        = isnan(s.Price);
Price(inan) = NaN;

% Returns with overnight
logret = diff(log(Price));
r      = ret(ret.Date == date & ret.Permno == permno,:);
if isnan(logret(1))
    logret(1) = r.RetCO;
else
    logret(1) = logret(1)+r.RetCO;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Den - variance of spyders
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Unsampled spy
unsampled         = getSpy(inf,date,date);
unsampled         = unsampled(~isInvalidTrade(unsampled),:);
[Datetime,~,subs] = unique(unsampled.Datetime);
Price             = accumarray(subs, unsampled.Price,[],@fast_median);

% Sample with MFE
Price = realized_price_filter(double(Price), Datetime,...
                                    'unit','fixed', grid+yyyymmdd2serial(date));

% Inherit same NaNs
Price(inan) = NaN;

% Returns with overnight
spylogret = diff(log(Price));
r         = ret(ret.Date == date & ret.Permno == 84398,:);
if isnan(spylogret(1))
    spylogret(1) = r.RetCO;
else
    spylogret(1) = spylogret(1)+r.RetCO;
end
num = nansum(logret.*spylogret);
den = nansum(spylogret.*spylogret);

record.ManualNum  = num;
record.ManualDen  = den;
record.ManualBeta = num./den;
disp(record)

%% Check if multiple days extends
fprintf('Checking beta with yearly lookback for a single.\n')

% Estimate
useovernight = false;
lookback     = 252;
betas        = getBetas(lookback,5,useovernight);

fprintf('%s: retrieving random series for manual check.\n', mfilename)

% Load spyders
spy = loadresults('spysampled5m');

% Load data
path2data = '.\data\TAQ\sampled\5min\';
if ~exist('master','var')
    master = load(fullfile(path2data,'master'),'-mat','mst');
end
idx  = ismember(master.mst.Permno,permno);
data = getTaqData(master.mst(idx,:), [],[],[],[],path2data);

fprintf('%s: calculating manual betas.\n', mfilename)

% Add spyders
data.Spy        = NaN(size(data,1),1);
[idx,pos]       = ismember(data.Datetime, spy.Datetime);
data.Spy(idx,1) = spy.Price(pos(idx));

% Get rid of NaNs
data = data(~isnan(data.Price) & ~isnan(data.Spy),:);

% Index by day
data.Date   = fix(data.Datetime);
manual      = table(unique(data.Date),'VariableNames',{'Date'});
N           = size(manual,1);
manual.Beta = NaN(N,1);

% Convert to return
data.Ret    = NaN(size(data,1),1);
data.Spyret = NaN(size(data,1),1);
for ii = 1:N
    idx              = data.Date == manual.Date(ii);
    prices           = data.Price(idx);
    spyprices        = data.Spy(idx);
    data.Ret(idx)    = [NaN; prices(2:end)./prices(1:end-1)-1];
    data.Spyret(idx) = [NaN; spyprices(2:end)./spyprices(1:end-1)-1];
end

% Drop overnight NaNs and empty days
data   = data(~isnan(data.Ret),:);
manual = manual(ismember(manual.Date, data.Date),:);

% Calculate Betas
N = size(manual,1);
for ii = lookback:N
    idx             = ismember(data.Date, manual.Date(ii-lookback+1:ii));
    num             = data.Ret(idx)'*data.Spyret(idx);
    den             = data.Spyret(idx)'*data.Spyret(idx);
    manual.Beta(ii) = num./den;
end
manual.Date = serial2yyyymmdd(manual.Date);

% Compare
compare = betas(betas.Permno == permno & ismember(betas.Date, manual.Date),:);

% Visual inspection
subplot(211)
title(sprintf('%d', permno))
plot(yyyymmdd2datetime(compare.Date), compare.Beta,...
     yyyymmdd2datetime(manual.Date), manual.Beta)
legend({'Automated','Manual'})

abdiff = abs(compare.Beta - manual.Beta);
counts = histc(abdiff, eps*10.^(0:15));
subplot(212)
bar((0:15)', counts)
set(gca,'XLim',[-1 16],'Xtick',0:15)
xlabel 'eps*10 to the nth power'