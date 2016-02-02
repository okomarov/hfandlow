%% Options
OPT_USEETF = false;

OPT_LAG    = 1;
OPT_PTF_UN = 5;

OPT_SHRINK = [0.4,0.6,0.6,0.6];
% OPT_SHRINK = [0.4,0.5,0.5,0.5];

%% Data
load('results\alldata_betaonly')

if OPT_USEETF
    etf           = reton(reton.Permno == 84398,:);
    etf           = etf(~isnan(etf.RetCC),:);
    [date,ia,ib]  = intersect(etf.Date, ff.Date);
    etf           = etf(ia,:);
    ff            = ff(ib,:);
    ff.MktMinusRF = etf.RetCC*100 - ff.RF;
    ret           = ret(ib,:);
    beta          = beta(ib,:,:);
end
%% Signals
[signals_LF, hpr, rf, mdate] = make_signals_LF(ret,date,ff,OPT_SHRINK);
signals_HF                   = make_signals_HF(xstr2num(permno),date,beta,OPT_SHRINK);

% NaN-out intersection
nsig = size(signals_LF,3);
inan = false(size(signals_LF));
for ii = 1:nsig
    inan(:,:,ii) = isnan(signals_LF(:,:,ii)) | isnan(signals_HF(:,:,ii));
end
signals_LF(inan) = NaN;
signals_HF(inan) = NaN;

%% Lag
if OPT_USEETF
    isMicro = isMicro(2:end,:);
end

% End-of-Month
signals_LF = signals_LF(1:end-OPT_LAG,:,:);
signals_HF = signals_HF(1:end-OPT_LAG,:,:);
isMicro    = isMicro(1:end-OPT_LAG,:);
% cap        = cap(1:end-OPT_LAG,:);

% Lag forward
hpr   = hpr(1+OPT_LAG:end,:);
rf    = rf(1+OPT_LAG:end,:);
mdate = mdate(1+OPT_LAG:end,:);

% Filter micro
hpr(isMicro) = NaN;
%% BAB ret
[ptfret, avgsig] = deal(cell(nsig,2));
for ii = 1:nsig
    [ptfret{ii,1},~,~,~,avgsig{ii,1}] = bab(hpr,signals_LF(:,:,ii),rf);
    [ptfret{ii,2},~,~,~,avgsig{ii,2}] = bab(hpr,signals_HF(:,:,ii),rf);
end
%% Desc stats
dt   = serial2datetime(datenum(1993,(1:size(signals_LF,1))+2,1)-1);
idec = dt >= yyyymmdd2datetime(20010501);

% Correlation
snames       = {'babm','babq','babs','baby','rbabm','rbabq','rbabs','rbaby'};
order        = reshape(reshape(1:nsig*2,nsig,2)',1,nsig*2);
allsig       = cat(3,signals_LF,signals_HF);
correlations = corrxs(allsig(:,:,order),snames(order));
correlations_dec = corrxs(allsig(idec,:,order),snames(order));

% Plot
figure
X    = datenum([min(dt(idec)), min(dt(idec)), max(dt(idec)),max(dt(idec))]);
YLIM = [0,6];
for ii = 1:nsig*2
    subplot(nsig,2,ii)
    inan      = isnan(ptfret{order(ii)});
    lvl       = cumprod(1+nan2zero(ptfret{order(ii)}));
    lvl(inan) = NaN;
    plot([min(dt(idec)),min(dt(idec))],YLIM,'-.','Color',[0.85,0.85,0.85],'LineWidth',1.5)
    hold on
    plot(dt,lvl)
    set(gca, 'TickLabelInterpreter','latex','Ylim',YLIM,'YTick',0:2:YLIM(2),'Layer','Top')
end

% Desc stats
desc   = cellfun(@(x) stratstats(dt,x*100,'Frequency','m','IsPercentageReturn',true), ptfret,'un',0);
desc   = cellfun(@(x) renameVarNames(x',{'Low','High','BAB'}), desc,'un',0);
catfun = @(sig,stats) [stats; array2table(nanmean(sig),'VariableNames',{'Low','High','BAB'},'RowNames',{'Avgsig'})];
desc   = cellfun(catfun,avgsig,desc,'un',0);
out = cell2mat(cellfun(@(x) x{:,:}, desc,'un',0));

% Desc stats, post decimalization
% Jan 2001 decimalization in NYSE and April in NASDAQ, take from May
desc2  = cellfun(@(x) stratstats(dt(idec),x(idec,:)*100,'Frequency','m','IsPercentageReturn',true), ptfret,'un',0);
desc2  = cellfun(@(x) renameVarNames(x',{'Low','High','BAB'}), desc2,'un',0);
catfun = @(sig,stats) [stats; array2table(nanmean(sig(idec,:)),'VariableNames',{'Low','High','BAB'},'RowNames',{'Avgsig'})];
desc2  = cellfun(catfun,avgsig,desc2,'un',0);
out = cell2mat(cellfun(@(x) x{:,:}, desc2,'un',0));
%% Tests
[coeff,se,tratio,pval] = regressHighOnLow(ptfret(:,2),ptfret(:,1));
[coeff2,se2,tratio2,pval2] = regressHighOnLow(ptfret(:,2),ptfret(:,1),idec);

[~,pValST,Zjk]   = cellfun(@(high,low) sharpetest(high(:,end), low(:,end)), ptfret(:,2),ptfret(:,1));
[~,pValST2,Zjk2] = cellfun(@(high,low) sharpetest(high(idec,end), low(idec,end)), ptfret(:,2),ptfret(:,1));
%% Risk-adjustment
factors = loadresults('RAfactors');
factors = factors(ismember(factors.Date, mdate),:);