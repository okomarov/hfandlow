%% BAB weights example
N    = 10;
z    = 1:N;
zbar = mean(z);
w    = [-min(z-zbar,0)/sum(abs(zbar-z))*2, max(z-zbar,0)/sum(abs(zbar-z))*2];

figure
set(gcf, 'Position', get(gcf,'Position').*[1,1,1,0.62],'PaperPositionMode','auto')
h = stem(z(:), reshape(w,N,[]), 'filled');
set(h(2),'Marker','diamond')
xticklabels = get(gca,'XTickLabel');
set(gca,'TickLabelInterpreter','latex','Xlim',[0,N+1],'XTickLabel',[' '; xticklabels;' '])
print('BabWeightExample','-depsc','-r200')
%% Strategy plots
ptfret = loadresults('ptfret_plot');

nsig  = size(ptfret,1);
order = reshape(reshape(1:nsig*2,nsig,2)',1,nsig*2);

dt      = serial2datetime(datenum(1993,(1:size(ptfret{1},1)+1)+1,1)-1);
idec    = dt >= yyyymmdd2datetime(20010501);
X       = datenum([min(dt(idec)), min(dt(idec)), max(dt(idec)),max(dt(idec))]);
YLIM    = [0,6];
mrkStep = 15;
for ii = 1:nsig*2
    figure
    set(gcf, 'Position', get(gcf,'Position').*[1,1,0.65,0.4],'PaperPositionMode','auto')
    inan      = [isnan(ptfret{order(ii)}); false(1,3)];
    lvl       = cumprod([zeros(1,3); nan2zero(ptfret{order(ii)})] +1);
    lvl(inan) = NaN;
    plot(dt,lvl)
    hold on
    h         = plot([min(dt(idec)),min(dt(idec))],YLIM,'-.','Color',[0.85,0.85,0.85],'LineWidth',1.5);
    uistack(h,'bottom')
    set(gca,'ColorOrderIndex',1)
    plot(dt(1:mrkStep:end),lvl(1:mrkStep:end,1),'x',dt(1:mrkStep:end),lvl(1:mrkStep:end,2),'+')
    set(gca, 'TickLabelInterpreter','latex','Ylim',YLIM,'YTick',0:2:YLIM(2),'Layer','Top')
    print(sprintf('strat%d',ii),'-depsc','-r200')
end
close all

YLIM  = [0,2];
decdt = [dt(find(idec,1,'first')-1), dt(idec)];
for ii = 1:nsig*2
    figure
    set(gcf, 'Position', get(gcf,'Position').*[1,1,0.65,0.4],'PaperPositionMode','auto')
    lvl = [ones(1,3); cumprod(1+ptfret{order(ii)}(idec(2:end),:))];
    plot(decdt,lvl)
    hold on
    h   = plot(get(gca,'Xlim'),[1,1],'-','Color',[0.5,0.5,0.5]);
    uistack(h,'bottom')
    set(gca,'ColorOrderIndex',1)
    plot(decdt(1:mrkStep:end),lvl(1:mrkStep:end,1),'x',decdt(1:mrkStep:end),lvl(1:mrkStep:end,2),'+')
    set(gca, 'TickLabelInterpreter','latex','Ylim',YLIM,'YTick',0:0.5:YLIM(2),'Layer','Top')
    print(sprintf('strat_dec%d',ii),'-depsc','-r200')
end
close all
%% Beta quintiles
% Justifying noise in beta estimates until decimalization
myunstack = @(tb,vname) sortrows(unstack(tb(:,{'Permno','Date',vname}),vname,'Permno'),'Date');

load('results\20160215_2237_betacomponents5m.mat')
res      = res(issp500member(res),:);
res.Beta = res.Num./res.Den;
res      = getMktCap(res,1);
res      = res(~isnan(res.Beta),:);
beta     = myunstack(res,'Beta');

% Monthly sampling
[~,pos] = unique(beta.Date/100,'last');
dates   = beta.Date(pos);
dt      = yyyymmdd2datetime(dates);
beta    = beta{pos,2:end};

% Beta quintiles
sz       = size(beta);
rows     = repmat((1:sz(1))',1,sz(2));
cols     = binSignal(beta);
betamean = accumarray([rows(:) cols(:)+1], beta(:), [],@nanmean);

% Quintiles sorted by cap
cap           = myunstack(res,'Cap');
cap           = cap{pos,2:end};
cols          = binSignal(cap);
betaMeanByCap = accumarray([rows(:) cols(:)+1], beta(:), [],@nanmean);

figure
set(gcf, 'Position', get(gcf,'Position').*[1,1,1,0.42],'PaperPositionMode','auto')
plot(dt, betamean(:,2:end))
h = recessionplot('recessions',datenum(reshape(dt([1,48,92,120]),2,2)'));
uistack(h,'bottom')
set(h,'FaceAlpha',1,'FaceColor',[0.9,0.9,0.9])
set(gca, 'TickLabelInterpreter','latex','Layer','Top','Ylim',[-1.5,2.5])
print('betaQuintiles','-depsc','-r200')

figure
set(gcf, 'Position', get(gcf,'Position').*[1,1,1,0.42],'PaperPositionMode','auto')
h = plot(dt, betaMeanByCap(:,[2,4,end]));
set(h,{'Color'},{[0    0.4470    0.7410]
               [0.9290    0.6940    0.1250]
               [0.4660    0.6740    0.1880]})
h = recessionplot('recessions',datenum(reshape(dt([1,48,92,120]),2,2)'));
uistack(h,'bottom')
set(h,'FaceAlpha',1,'FaceColor',[0.9,0.9,0.9])
set(gca, 'TickLabelInterpreter','latex','Layer','Top','Ylim',[-0.4,1.6])
print('betaQuintilesCap','-depsc','-r200')
%% Count 0 returns

testname = 'countnullrets';
try
    counts = loadresults(testname);
catch
    path2data = '..\data\TAQ\sampled\5min\nobad_vw';
    counts    = AnalyzeHflow(testname,[],[], path2data,[],8);
end

counts = counts(~isprobdate(counts.Date),:);

% Monthly sampling
[~,pos] = unique(counts.Date/100,'last');
dates   = counts.Date(pos);
dt      = yyyymmdd2datetime(dates);

% Only sp500 members
isp500    = issp500member(counts);
[~,subs]  = ismember(counts.Date(isp500)/100,dates/100);
avgcounts = accumarray(subs, counts.Nullrets(isp500),[],@mean);

figure
set(gcf, 'Position', get(gcf,'Position').*[1,1,1,0.42],'PaperPositionMode','auto')
plot(dt, avgcounts)
h = recessionplot('recessions',datenum(reshape(dt([1,48,92,120]),2,2)'));
uistack(h,'bottom')
set(h,'FaceAlpha',1,'FaceColor',[0.9,0.9,0.9])
set(gca, 'TickLabelInterpreter','latex','Layer','Top','Ylim',[0,80])
print('countNullConst','-depsc','-r200')

% Spyders
spy = sortrows(counts(counts.Permno == 84398,{'Date','Nullrets'}),'Date');

[~,subs]  = ismember(spy.Date/100,dates/100);
avgcounts = accumarray(subs, spy.Nullrets,[],@mean);
avgcounts(1) = NaN;

figure
set(gcf, 'Position', get(gcf,'Position').*[1,1,1,0.42],'PaperPositionMode','auto')
plot(dt, avgcounts)
h = recessionplot('recessions',datenum(reshape(dt([1,48,92,120]),2,2)'));
uistack(h,'bottom')
set(h,'FaceAlpha',1,'FaceColor',[0.9,0.9,0.9])
set(gca, 'TickLabelInterpreter','latex','Layer','Top','Ylim',[0,80])
print('countNullSpy','-depsc','-r200')

%% Trends in SP500 betas
myunstack = @(tb,vname) sortrows(unstack(tb(:,{'Permno','Date',vname}),vname,'Permno'),'Date');

load('results\20160114_1732_betacomponents5mon.mat')
res      = res(issp500member(res),:);
res.Beta = res.Num./res.Den;
res      = res(~isnan(res.Beta),:);
beta     = myunstack(res,'Beta');
plot(prctile(beta{1:20:end,2:end},10:10:90,2))

load('results\20160118_1949_betacomponents5m.mat')
res      = res(issp500member(res),:);
res.Beta = res.Num./res.Den;
res      = res(~isnan(res.Beta),:);
beta     = myunstack(res,'Beta');
plot(prctile(beta{1:20:end,2:end},10:10:90,2))

% You can't see the effect of decimalization on beta75 minutes
load('results\20160121_1655_betacomponents75mon.mat')
res      = res(issp500member(res),:);
res.Beta = res.Num./res.Den;
res      = res(~isinf(res.Beta),:);
res      = res(~isnan(res.Beta),:);
beta     = myunstack(res,'Beta');
plot(prctile(beta{1:20:end,2:end},10:10:90,2))

mean(nanmean(beta{1:2046,2:end},2))
mean(nanmedian(beta{1:2046,2:end},2))
mean(nanmean(beta{2046:end,2:end},2))
mean(nanmedian(beta{2046:end,2:end},2))

load('results\20160121_2210_betacomponents75m.mat')
res      = res(issp500member(res),:);
res.Beta = res.Num./res.Den;
res      = res(~isinf(res.Beta),:);
res      = res(~isnan(res.Beta),:);
beta     = myunstack(res,'Beta');
plot(prctile(beta{1:20:end,2:end},10:10:90,2))
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
