%% Options
OPT_VW = true;

OPT_LAG    = 1;
OPT_PTF_UN = 5;

OPT_SHRINK = [0.4,0.5,0.5,0.5];

%% Data
load('results\alldata_betaonly')

%% Signals

% Low freqeuncy signals
[signals_LF, hpr, rf, mdate] = make_signals_LF(ret,date,ff,OPT_SHRINK);

% High freqeuncy signals
signals_HF = make_signals_HF(xstr2num(permno),date,beta,OPT_SHRINK);
nsig       = size(signals_LF,3);
inan       = false(size(signals_LF));
for ii = 1:nsig
    inan(:,:,ii) = isnan(signals_LF(:,:,ii)) | isnan(signals_HF(:,:,ii));
end
signals_LF(inan) = NaN;
signals_HF(inan) = NaN;

snames       = {'babm','babq','babs','baby','rbabm','rbabq','rbabs','rbaby'};
order        = reshape(reshape(1:nsig*2,nsig,2)',1,nsig*2);
allsig       = cat(3,signals_LF,signals_HF);
correlations = corrxs(allsig(:,:,order),snames(order));
%% Lag
% End-of-Month
signals_LF = signals_LF(1:end-OPT_LAG,:,:);
signals_HF = signals_HF(1:end-OPT_LAG,:,:);
isMicro    = isMicro(1:end-OPT_LAG,:);
% cap        = cap(1:end-OPT_LAG,:);

% Lag forward
hpr   = hpr(1+OPT_LAG:end,:);
rf    = rf(1+OPT_LAG:end,:);
mdate = mdate(1+OPT_LAG:end,:);
%% Filter micro
hpr(isMicro) = NaN;
%% PTFRET

[ptfret, avgsig] = deal(cell(nsig,2));
% Bab
for ii = 1:nsig
    [ptfret{ii,1},~,~,~,avgsig{ii,1}] = bab(hpr,signals_LF(:,:,ii),rf);
    [ptfret{ii,2},~,~,~,avgsig{ii,2}] = bab(hpr,signals_HF(:,:,ii),rf);
end

dt     = serial2datetime(datenum(1993,(1:size(hpr,1))+2,1)-1);
desc   = cellfun(@(x) stratstats(dt,x*100,'Frequency','m','IsPercentageReturn',true), ptfret,'un',0);
desc   = cellfun(@(x) renameVarNames(x',{'Low','High','BAB'}), desc,'un',0);
catfun = @(sig,stats) [stats; array2table(nanmean(sig),'VariableNames',{'Low','High','BAB'},'RowNames',{'Avgsig'})];
desc   = cellfun(catfun,avgsig,desc,'un',0);

figure
for ii = 1:nsig*2
    subplot(nsig,2,ii)
    inan      = isnan(ptfret{order(ii)});
    lvl       = cumprod(1+nan2zero(ptfret{order(ii)}));
    lvl(inan) = NaN;
    plot(dt,lvl)
    set(gca, 'TickLabelInterpreter','latex','Ylim',[0,8],'YTick',0:2:8)
end

%% PTFRET after decimalization
% Jan 2001 decimalization in NYSE and April in NASDAQ, take from May

dt   = serial2datetime(datenum(1993,(1:size(hpr,1))+2,1)-1);
idec = dt >= yyyymmdd2datetime(20010501);

desc2  = cellfun(@(x) stratstats(dt(idec),x(idec,:)*100,'Frequency','m','IsPercentageReturn',true), ptfret,'un',0);
desc2  = cellfun(@(x) renameVarNames(x',{'Low','High','BAB'}), desc2,'un',0);
catfun = @(sig,stats) [stats; array2table(nanmean(sig(idec,:)),'VariableNames',{'Low','High','BAB'},'RowNames',{'Avgsig'})];
desc2  = cellfun(catfun,avgsig,desc2,'un',0);


%% Tests
[coeff,se] = deal(NaN(nsig,2));
for ii = 1:nsig
    inan = isnan(ptfret{ii,1});
    nobs = nnz(~inan);
    opts = {'bandwidth',floor(4*(nobs/100)^(2/9))+1,'intercept',false,'type','HAC','weights','BT','display','off'};
    f    = @(x,y) hac(x, y, opts{:});

    l                        = ones(size(ptfret{ii,2},1),1);
    [~,se(ii,:),coeff(ii,:)] = f([l, ptfret{ii,1}(:,end)], ptfret{ii,2}(:,end));
end
% tstat and pvalues
tratio = coeff./se;
pval   = 2 * normcdf(-abs(tratio));

[coeff,se] = deal(NaN(nsig,2));
for ii = 1:nsig
    nobs = nnz(idec);
    opts = {'bandwidth',floor(4*(nobs/100)^(2/9))+1,'intercept',false,'type','HAC','weights','BT','display','off'};
    f    = @(x,y) hac(x, y, opts{:});

    l                        = ones(size(ptfret{ii,2}(idec,:),1),1);
    [~,se(ii,:),coeff(ii,:)] = f([l, ptfret{ii,1}(idec,end)], ptfret{ii,2}(idec,end));
end
% tstat and pvalues
tratio = coeff./se;
pval   = 2 * normcdf(-abs(tratio));

sharpetest(ptfret{1,2}(idec,end), ptfret{1,1}(idec,1))
%% Risk-adjustment
factors = loadresults('RAfactors');
factors = factors(ismember(factors.Date, mdate),:);
%% RV, RCOV, RBETA plots

% Alcoa
aa          = getTaqData('symbol','AA',[],[],[],'..\data\TAQ\sampled\5min\nobad_vw');
% issorted(aa.Datetime)
aa.Ret      = [NaN; diff(log(aa.Price))];
ion         = [true; diff(fix(aa.Datetime))~=0];
aa.Ret(ion) = NaN;

spy          = getSpy(5);
spy.Ret      = [NaN; diff(log(spy.Price))];
ion          = [true; diff(fix(spy.Datetime))~=0];
spy.Ret(ion) = NaN;

% Intersect dates
Date = intersect(fix(aa.Datetime), fix(spy.Datetime));
aa   = aa(ismember(fix(aa.Datetime), Date),:);
spy  = spy(ismember(fix(spy.Datetime), Date),:);

% NO OVERNIGHT
%%%%%%%%%%%%%%
rv   = @(subs) sqrt(accumarray(subs,spy.Ret.^2,[],@nansum));
rcov = @(subs,stock) accumarray(subs,spy.Ret.*stock.Ret,[],@nansum);

% Daily
[Date,~,iday] = unique(serial2yyyymmdd(spy.Datetime));
RVd           = rv(iday);
RCVd          = rcov(iday,aa);

% Monthly
[~,pos,imonth] = unique(serial2yyyymmdd(spy.Datetime)/100,'last');
RVm            = rv(imonth);
RCVm           = rcov(imonth,aa);

% Yearly, rolling every month
RVy      = NaN(size(RCVm));
RCVy     = NaN(size(RCVm));
rcov_run = @(idx,stock) nansum(spy.Ret(idx).*stock.Ret(idx));
for ii = 12:size(RCVm,1)
    idx      = ismember(imonth, ii-11:ii);
    RVy(ii)  = sqrt(nansum(spy.Ret(idx).^2));
    RCVy(ii) = rcov_run(idx,aa);
end

% WITH OVERNIGHT
%%%%%%%%%%%%%%%%
spy.Permno(:,1) = 84398;
spy             = addOvernightRet(spy);
aa              = addOvernightRet(aa);

rv   = @(subs) sqrt(accumarray(subs,spy.Ret.^2,[],@nansum));
rcov = @(subs,stock) accumarray(subs,spy.Ret.*stock.Ret,[],@nansum);

RVd_on  = rv(iday);
RCVd_on = rcov(iday,aa);

RVm_on  = rv(imonth);
RCVm_on = rcov(imonth,aa);

RVy_on   = NaN(size(RVm));
RCVy_on  = NaN(size(RVm));
rcov_run = @(idx,stock) nansum(spy.Ret(idx).*stock.Ret(idx));
for ii = 12:size(RVm,1)
    idx         = ismember(imonth, ii-11:ii);
    RVy_on(ii)  = sqrt(nansum(spy.Ret(idx).^2));
    RCVy_on(ii) = rcov_run(idx,aa);
end

dt_daily   = yyyymmdd2datetime(Date);
dt_monthly = serial2datetime(spy.Datetime(pos));

% Plot RV
figure
set(gcf, 'Position', get(gcf,'Position').*[1,1,1,0.62],'PaperPositionMode','auto')
plot(dt_daily, RVd,    dt_monthly, RVm,    dt_monthly, RVy,...
     dt_daily, RVd_on, dt_monthly, RVm_on, dt_monthly, RVy_on)
legend('daily','monthly','yearly','d_on','m_on','y_on','Location','NorthWest')
legend boxoff
set(gca,'Ytick',[0,0.1,0.2,0.3,0.4],'YtickLabel',{'0','10\%','20\%','30\%','40\%'})
set(gca,'TickLabelInterpreter','latex')
print('SPYrv','-depsc','-r200')

% Plot RCV
figure
set(gcf, 'Position', get(gcf,'Position').*[1,1,1,0.62],'PaperPositionMode','auto')
plot(dt_daily, RCVd,    dt_monthly, RCVm,    dt_monthly, RCVy,...
     dt_daily, RCVd_on, dt_monthly, RCVm_on, dt_monthly, RCVy_on)
legend('daily','monthly','yearly','d_on','m_on','y_on','Location','NorthWest')
legend boxoff
set(gca,'Ytick',[0,0.1,0.2,0.3],'YtickLabel',{'0','10\%','20\%','30\%'},'Ylim',[-0.01,0.3])
set(gca,'TickLabelInterpreter','latex')

% Plot RBETA
figure
set(gcf, 'Position', get(gcf,'Position').*[1,1,1,0.62],'PaperPositionMode','auto')
plot(dt_daily, RCVd./RVd.^2,    dt_monthly, RCVm./RVm.^2,    dt_monthly, RCVy./RVy.^2,...
     dt_daily, RCVd_on./RVd_on.^2, dt_monthly, RCVm_on./RVm_on.^2, dt_monthly, RCVy_on./RVy_on.^2)
legend('daily','monthly','yearly','d_on','m_on','y_on','Location','NorthWest')
legend boxoff
set(gca,'Ytick',[0,0.1,0.2,0.3],'YtickLabel',{'0','10\%','20\%','30\%'},'Ylim',[-0.01,0.3])
set(gca,'TickLabelInterpreter','latex')
