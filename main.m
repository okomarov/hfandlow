%% Options
OPT_VW = true;

OPT_LAG    = 1;
OPT_PTF_UN = 5;

%% Data
load('results\alldata')
master = loadresults('master');
reton  = loadresults('reton');

% % Things in alldata, already unstacked
% dsf         = loadresults('dsf');
% dsf.IsMicro = isMicrocap(dsf,'Prc');
% dsf         = getMktCap(dsf);
% rskew       = loadresults('skew');
% beta        = loadresults('beta5minon');
% ff          = loadresults('F-F_Research_Data_5_Factors_2x3_daily_TXT');
%% Signals

% Low freqeuncy signals
[signals_LF, hpr, rf, mdate] = make_signals_LF(ret,date,ff);

% High freqeuncy signals
signals_HF = make_signals_HF(xstr2num(permno),date,master,reton,ff,rskew,beta);

inanLF = isnan(signals_LF);
inanHF = isnan(signals_HF);
inan = cat(3,repmat(any(inanHF(:,:,1:3) | inanLF(:,:,1:3),3),[1,1,3]),...
                    any(inanHF(:,:,4)   | inanLF(:,:,4)  ,3));
% plot(sum(~any(inanHF(:,:,4)   | inanLF(:,:,4)  ,3),2))  
% plot(sum(~any(inanHF(:,:,1:3) | inanLF(:,:,1:3),3),2))
signals_HF(inan) = NaN;
signals_LF(inan) = NaN;

snames = {'ca','rskd','hsk','bab','rca','rskd5','rskm5','rbab'};
order  = [1,5,3,2,7,6,4,8];
allsig = cat(3,signals_LF,signals_HF);
correlations = corrxs(allsig(:,:,order),snames(order));
%% Lag
% End-of-Month
signals_LF = signals_LF(1:end-OPT_LAG,:,:);
signals_HF = signals_HF(1:end-OPT_LAG,:,:);
isMicro    = isMicro(1:end-OPT_LAG,:);
cap        = cap(1:end-OPT_LAG,:);

% Lag forward
hpr   = hpr(1+OPT_LAG:end,:);
rf    = rf(1+OPT_LAG:end,:);
mdate = mdate(1+OPT_LAG:end,:);
%% Filter micro
hpr(isMicro) = NaN;
%% PTFRET
if OPT_VW
    opts = struct('PortfolioNumber',OPT_PTF_UN, 'Weights',double(cap));
else
    opts = struct('PortfolioNumber',OPT_PTF_UN);
end

% Alpha
[ptfret{1,1},~,counts{1,1},avgsig{1,1}] = portfolio_sort(hpr, signals_LF(:,:,1), opts);
[ptfret{1,2},~,counts{1,2},avgsig{1,2}] = portfolio_sort(hpr, signals_HF(:,:,1), opts);

% Skewness
[ptfret{2,1},~,counts{2,1},avgsig{2,1}] = portfolio_sort(hpr, signals_LF(:,:,2), opts);
[ptfret{2,2},~,counts{2,2},avgsig{2,2}] = portfolio_sort(hpr, signals_HF(:,:,2), opts);

% Skewness#2
[ptfret{3,1},~,counts{3,1},avgsig{3,1}] = portfolio_sort(hpr, signals_LF(:,:,3), opts);
[ptfret{3,2},~,counts{3,2},avgsig{3,2}] = portfolio_sort(hpr, signals_HF(:,:,3), opts);

% Bab
[ptfret{4,1},~,~,~,avgsig{4,1}] = bab(hpr,signals_LF(:,:,4),rf);
[ptfret{4,2},~,~,~,avgsig{4,2}] = bab(hpr,signals_HF(:,:,4),rf);

dt = serial2datetime(datenum(1993,(1:size(hpr,1))+2,1)-1);
desc = cellfun(@(r) stratstats(dt,r*100,'Frequency','m','IsPercentageReturn',true), ptfret,'un',0);

figure
for r = 1:4
    for c = 1:2
        n = (r-1)*2+c;
        subplot(420+n)
        if r < 4
            plot(dt,cumprod(1+ptfret{r,c}))
        else
            plot(dt(12:end),cumprod(1+ptfret{r,c}(12:end,:)))
        end
        axis tight
        set(gca, 'Ylim',[0,10])
    end
end
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
