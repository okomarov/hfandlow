
OPT_LAGDAY   = 1;
OPT_BETAFREQ = 75;
OPT_USEON    = true;

switch OPT_BETAFREQ
    case 75
        grid = [0, 9.75/24:OPT_BETAFREQ/(60*24):16/24];
    case 30
        grid = [0, 10/24:OPT_BETAFREQ/(60*24):16/24];
    otherwise 
        grid = [];
end
%% Select data
% Index data
datapath = '..\data\TAQ\sampled\5min\nobad_vw';
master   = load(fullfile(datapath,'master'),'-mat');
master   = master.mst(master.mst.Permno ~= 0,:);
master   = sortrows(master,{'Permno','Date'});

% Get market
mkt = master(master.Permno == 84398,:);

% Common shares
idx    = iscommonshare(master);
master = master(idx,:);

% Incomplete days
idx    = isprobdate(master.Date);
master = master(~idx,:);

% Minobs
res               = loadresults('countBadPrices','..\results');
isEnoughObs       = (res.Ntot - res.Nbadtot) >= 79;
res               = addPermno(res);
[~,pos]           = ismembIdDate(master.Permno, master.Date, res.Permno, res.Date);
isEnoughObs       = isEnoughObs(pos,:);
[~,idx,pos]       = lagpanel(master(:,{'Date','Permno'}),'Permno');
isEnoughObs(~idx) = isEnoughObs(pos);
master            = master(isEnoughObs,:);

% % Count
% [~,~,subs] = unique(master.Date);
% accumarray(subs,1)

% Beta components at 75, i.e. return from open to 11:00, 12:15, 13:30, 14:45, 16:00
try
    if OPT_USEON
        beta = loadresults(sprintf('betacomponents%dmon',OPT_BETAFREQ));
    else
        beta = loadresults(sprintf('betacomponents%dm',OPT_BETAFREQ));
    end
catch
    beta = estimateBetaComponents(OPT_BETAFREQ,OPT_USEON,false,grid);
end
% Drop january, no beta spyders, hence no betacomponents
beta = beta(beta.Date/100 ~= 199301, :);

[~,ia,ib] = intersectIdDate(beta.Permno, beta.Date,master.Permno, master.Date);
beta      = beta(ia,:);
master    = master(ib,:);

% CRSP returns
dsf       = loadresults('dsfquery','..\results');
[~,ia,ib] = intersectIdDate(dsf.Permno, dsf.Date,master.Permno, master.Date);
dsf       = dsf(ia,:);
master    = master(ib,:);

% Beta components - re-run
idx  = ismembIdDate(beta.Permno, beta.Date, master.Permno, master.Date);
beta = beta(idx,:);

% Add back mkt
master = [master; mkt];

% % Overnight returns
% reton = loadresults('return_intraday_overnight');
% idx   = ismembIdDate(reton.Permno, reton.Date, master.Permno, master.Date);
% reton = reton(idx,:);

% importFrenchData('F-F_Research_Data_5_Factors_2x3_daily_TXT.zip','results');
%% Second stage
myunstack = @(tb,vname) sortrows(unstack(tb(:,{'Permno','Date',vname}),vname,'Permno'),'Date');

% Returns
ret    = myunstack(dsf,'Ret');
date   = ret.Date;
permno = ret.Properties.VariableNames(2:end);
ret    = double(ret{:,2:end});

% End-of-Month ismicro
dsf.IsMicro = isMicrocap(dsf,'Prc');
isMicro     = myunstack(dsf,'IsMicro');
[~,pos]     = unique(isMicro.Date/100,'last');
isMicro     = isMicro{pos,2:end};

% Beta components
num               = myunstack(beta,'Num');
den               = myunstack(beta,'Den');
beta              = cat(3,num{:,2:end},den{:,2:end});
beta(isinf(beta)) = NaN; % permno 46288 on 19931008 is delisted with close-to-close return of -100%
clear den num

% Factors
ff = loadresults('F-F_Research_Data_5_Factors_2x3_daily_TXT');
ff = ff(ismember(ff.Date, unique(dsf.Date)),:);

save(sprintf('results\\alldata_beta%d',OPT_BETAFREQ), 'master', 'date', 'permno', 'ret', 'isMicro', 'beta', 'ff')
%% Illiquidity
myunstack = @(tb,vname) sortrows(unstack(tb(:,{'Permno','Date',vname}),vname,'Permno'),'Date');
load('results\alldata_beta75','permno','date')

dsf = loadresults('dsfquery','..\results');
idx = ismember(dsf.Permno, xstr2num(permno));
dsf = dsf(idx,:);

% Incomplete days
idx = isprobdate(dsf.Date);
dsf = dsf(~idx,:);

% Amihud illiq |r_t|/(price_t.*volume_t)
dsf.Amihud       = abs(dsf.Ret ./ (double(dsf.Prc).*double(dsf.Vol/100)));
iinf             = isinf(dsf.Amihud);
dsf.Amihud(iinf) = NaN;
illiq            = myunstack(dsf,'Amihud');

idx   = ismember(illiq.Date, date);
illiq = illiq(idx,:);

% End-of-month rolling average of a year of daily observations
[mdate,pos,midx] = unique(illiq.Date/100,'last');
illiq            = illiq{:,2:end};

nmonths = numel(mdate);
out     = NaN(nmonths,numel(permno));
for ii = 12:nmonths
    idx       = ismember(midx, ii-12+1:ii);
    out(ii,:) = nanmean(illiq(idx,:),1);
end
illiq = log(out);
% illiq = out;
save results\illiq illiq