
OPT_LAGDAY   = 1;
OPT_BETAFREQ = 75;
OPT_USEON = false;
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
res         = loadresults('countBadPrices','..\results');
isEnoughObs = (res.Ntot - res.Nbadtot) >= 79;
[~,pos]     = ismembIdDate(master.Id, master.Date, res.Id, res.Date);
isEnoughObs = isEnoughObs(pos,:);
isEnoughObs = [false(OPT_LAGDAY,1); isEnoughObs(1:end-OPT_LAGDAY)];
master      = master(isEnoughObs,:);

% % Count
% tmp       = sortrows(unstack(master(:,{'Date','Permno','File'}), 'File','Permno'),'Date');
% tmp       = tmp{:,2:end};
% count_all = sum(tmp~=0,2);

% Beta components at 75, i.e. return from open to 11:00, 12:15, 13:30, 14:45, 16:00
try
    beta = loadresults(sprintf('betacomponents%dmon',OPT_BETAFREQ));
catch
    grid        = [0 11/24:OPT_BETAFREQ/(60*24):16/24];
    half_second = 0.5/(60*60*24);
    grid        = grid + half_second;
    beta        = estimateBetaComponents(OPT_BETAFREQ,true,false,grid);
end
[~,ia,ib] = intersectIdDate(beta.Permno, beta.Date,master.Permno, master.Date);
beta      = beta(ia,:);
master    = master(ib,:);

% CRSP returns
dsf       = loadresults('dsfquery','..\results');
[~,ia,ib] = intersectIdDate(dsf.Permno, dsf.Date,master.Permno, master.Date);
dsf       = dsf(ia,:);
master    = master(ib,:);

% Skewness
try
    skew = loadresults('skewcomponents');
catch
    reton                    = loadresults('return_intraday_overnight');
    [idx,pos]                = ismembIdDate(reton.Permno, reton.Date,master.Permno, master.Date);
    master.RetCO(pos(idx),1) = reton.RetCO(idx);
    mst                      = cache2cell(master,master.File);
    skew                     = AnalyzeHflow('skewcomponents',[],mst,datapath,[],8);
end
[idx,~] = ismembIdDate(skew.Permno, skew.Date, master.Permno, master.Date);
skew    = skew(idx,:);

% Beta components - re-run
idx  = ismembIdDate(beta.Permno, beta.Date, master.Permno, master.Date);
beta = beta(idx,:);

% Add back mkt
master = [master; mkt];

% Overnight returns
reton = loadresults('return_intraday_overnight');
idx   = ismembIdDate(reton.Permno, reton.Date, master.Permno, master.Date);
reton = reton(idx,:);

save('results\dsf.mat','dsf')
save('results\master.mat','master')
save(fullfile('results',sprintf('beta%dminon.mat',OPT_BETAFREQ)),'beta')
save('results\skew.mat','skew')
save('results\reton.mat','reton')

% importFrenchData('F-F_Research_Data_5_Factors_2x3_daily_TXT.zip','results');
%% Second stage
myunstack = @(tb,vname) sortrows(unstack(tb(:,{'Permno','Date',vname}),vname,'Permno'),'Date');

% Returns
dsf    = loadresults('dsf');
ret    = myunstack(dsf,'Ret');
date   = ret.Date;
permno = ret.Properties.VariableNames(2:end);
ret    = double(ret{:,2:end});

% End-of-Month ismicro
dsf.IsMicro = isMicrocap(dsf,'Prc');
isMicro     = myunstack(dsf,'IsMicro');
[~,pos]     = unique(isMicro.Date/100,'last');
isMicro     = isMicro{pos,2:end};

% Mkt cap
dsf     = getMktCap(dsf);
cap     = myunstack(dsf,'Cap');
[~,pos] = unique(cap.Date/100,'last');
cap     = cap{pos,2:end};

% Realized skewness
rskew   = loadresults('skew');
tmp.N   = myunstack(rskew,'N');
tmp.Sx3 = myunstack(rskew,'Sx3');
tmp.Rv  = myunstack(rskew,'Rv');
tmp.N   = tmp.N{:,2:end};
tmp.Sx3 = tmp.Sx3{:,2:end};
tmp.Rv  = tmp.Rv{:,2:end};
rskew   = tmp;
clear tmp

% Overnight return
reton = loadresults('reton');

% Beta components
beta              = loadresults(sprintf('beta%dminon',OPT_BETAFREQ));
num               = myunstack(beta,'Num');
den               = myunstack(beta,'Den');
beta              = cat(3,num{:,2:end},den{:,2:end});
beta(isinf(beta)) = NaN; % permno 46288 on 19931008 is delisted with close-to-close return of -100%
clear den num

% Factors
ff = loadresults('F-F_Research_Data_5_Factors_2x3_daily_TXT');
ff = ff(ismember(ff.Date, unique(dsf.Date)),:);

save results\alldata date permno ret isMicro cap rskew beta ff