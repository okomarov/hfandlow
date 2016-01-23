function sp500proxy = sp500intraday(path2data,freq)
% ISSUES:
% - Should use previous day close of the index, and update it as soon as
% constituents start to trade. However, it requires split/merge/divid
% adjustment. Now, using back-filling.
% - Filtering out by membership should be complemented with dragging of the
% last avilable number of shares


cap = loadresults('mktcap');
cap = cap(issp500member(cap),:);

taq2crsp  = loadresults('taq2crsp_sp500');
spconst   = loadresults('spconst');
dseshares = loadresults('dseshares');
if nargin < 1 || isempty(path2data), path2data = '.\data\TAQ\sampled\5min'; end
master = load(fullfile(path2data, 'master'), '-mat');

% Filter out non members 
dseshares  = dseshares(ismember(dseshares.PERMNO, spconst.Id),:);
isp500     = issp500member(master.mst(:,{'Date','UnID'}));
master.mst = master.mst(isp500,:);

% Time consolidation
idx       = isfeatchange(dseshares(:,{'PERMNO', 'SHROUT','SHRSDT'}));
from      = find(idx);
to        = [from(2:end)-1; numel(idx)];
dseshares = [dseshares(from,{'PERMNO','SHROUT','SHRSDT'}), dseshares(to,'SHRENDDT')];
dseshares = pivotFromTo(dseshares(:,{'PERMNO','SHRSDT','SHRENDDT','SHROUT'}));

% Sample at reference dates
refdates        = unique(master.mst.Date);
dseshares.Panel = sampledates(dseshares.Panel,refdates);
spconst.Panel   = sampledates(spconst.Panel,refdates);

% Filter out by membership
vnames = getVariableNames(spconst.Panel);
for c = 2:size(spconst.Panel,2)
    id = vnames{c};
    dseshares.Panel.(id)(~spconst.Panel.(id)) = 0;
end

% Data matrix of number of shares
sharesdata = table2array(dseshares.Panel);

% Add permno to master
[~,pos]           = ismember(master.mst.UnID,taq2crsp.ID);
master.mst.Permno = taq2crsp.permno(pos);

% LOOP by date
[dates,~,subs] = unique(master.mst.Date);
ndates         = numel(dates);
sp500proxy     = cell(ndates,1);
tic
poolStartup(4, 'AttachedFiles',{'.\utils\poolStartup.m'})
pctRunOnAll warning off MATLAB:table:ModifiedVarnames
parfor ii = 1:ndates
    disp(ii)
    % Retrieve reference permnos for given date
    row    = ismembc2(dates(ii), spconst.Panel.Date);
    shares = sharesdata(row,2:end);
    
    % Get mst records
    idates = subs == ii;
    mst    = master.mst(idates,:);
    
    % Restrict to reference permnos (overlapping?)
    [~, imst, ishares] = intersect(mst.Permno, spconst.Id);
    mst = mst(imst,:);
    
    % Get data
    data = getTaqData(mst,[],[],[],[],path2data);
    
    % Use the backfill price strategy to build proxy (should update
    % previous day index with gradually changing prices, need overnight return)
    data.Price = flipud(nanfillts(data.Price(end:-1:1)));
    
    % Ensure all datetimes are to the minute
    dt = datevec(data.Datetime);
    data.Datetime = datenum(dt(:,1),dt(:,2),dt(:,3),dt(:,4),dt(:,5),0);
    
    % Pivot
    prices     = unstack(data(:,{'Permno','Datetime','Price'}),'Price','Permno');
    datetimes  = prices.Datetime;
    prices     = table2array(prices(:,2:end));
    
    % (should use pervious day close)
    capitaliz = double(shares(ishares)).* prices(1,:);
    index = sum(bsxfun(@times, prices, capitaliz./sum(capitaliz)),2);
    
    % Store results
    sp500proxy{ii} = table(datetimes, index, 'VariableNames',{'Datetime','Price'});
end
pctRunOnAll warning on MATLAB:table:ModifiedVarnames
delete(gcp)
sp500proxy = cat(1,sp500proxy{:});

% Save
matname = sprintf('%s_sp500proxy%dm.mat',datestr(now,'yyyymmdd_HHMM'),freq);
save(fullfile('.\results',matname), 'sp500proxy')

%% Plot vs spyders
if nargout == 0
    spysampled = loadresults('spysampled');
    
    subplot(211)
    plot(datetime(datevec(sp500proxy.Datetime)), sp500proxy.Price)
    title 'sp500 proxy - rebalanced daily with open price'
    
    subplot(212)
    plot(datetime(datevec(spysampled.Datetime)), spysampled.Price)
    title 'spyders'
    
    saveas(gcf,'.\results\SP500proxy vs SPY.png')
end
toc
end