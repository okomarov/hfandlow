function betaPercentiles(ptiles, lookback, freq, useon, useproxy, issp, iscs)
% betaPercentiles(ptiles, lookback, freq, useon, useproxy, issp, iscs)

if nargin < 1 || isempty(ptiles),   ptiles   = 10:10:90;end
if nargin < 2 || isempty(lookback), lookback = 1;       end
if nargin < 3 || isempty(freq),     freq     = 5;       end
if nargin < 4 || isempty(useon),    useon    = true;    end
if nargin < 5 || isempty(useproxy), useproxy = false;   end
if nargin < 6 || isempty(issp),     issp     = true;   end
if nargin < 7 || isempty(iscs),     iscs     = true;   end

% Get Betas
betas = getBetas(lookback, freq, useon, useproxy, issp, iscs);

% Filter out problematic dates
betas = betas(~isprobdate(betas.Date),:);

% Sample/expand
refdates = serial2yyyymmdd(datenum(1993,2:234,1)-1);
betas    = sampledates(betas,refdates,1);

% All days
% refdates = Betas.Date;
% tmp      = table2array(Betas);

% Plot
plotdates   = yyyymmdd2datetime(refdates);
percentiles = prctile(table2array(betas(:,2:end)),ptiles,2);
plot(plotdates, percentiles)
legend(arrayfun(@(x) sprintf('%d^{th} ',x),ptiles,'un',0))

name = matname('BetaPercentiles');
title({'Cross-sectional percentiles of';  name},'interpreter','none');
saveas(gcf, fullfile('results','fig',sprintf('%s.png',name)))

    function name = matname(name)
        if useon,    useon    = 'withON';        else useon    = 'noON'; end
        if useproxy, useproxy = 'proxy';         else useproxy = 'spy'; end
        if issp,     issp     = 'sp500';         else issp     = 'allTAQ'; end
        if iscs,     iscs     = 'commonshares';  else iscs     = 'allshares'; end
        name = sprintf('%s_%dm%dd_%s_%s_%s_%s', name, freq, lookback,useon, useproxy,issp,iscs);
    end

end
