%% Options
OPT_USEETF = false;

OPT_FREQ = 75;

OPT_LAG    = 1;
OPT_PTF_UN = 5;
OPT_PTF_DB = 5;

OPT_SHRINK = [0.4,0.6,0.6,0.6];
% OPT_SHRINK = [0.8,0.8,0.8,0.8];

OPT_BLOCKS_DEC = {1 6 2 1}'; % whole horizon
% OPT_BLOCKS_DEC = {1 8 2 1}';
%% Data
load(sprintf('results\\alldata_beta%d',OPT_FREQ))

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
nsig     = size(signals_LF,3);
inan_sig = false(size(signals_LF));
for ii = 1:nsig
    inan_sig(:,:,ii) = isnan(signals_LF(:,:,ii)) | isnan(signals_HF(:,:,ii));
end
signals_LF(inan_sig) = NaN;
signals_HF(inan_sig) = NaN;

%% Lag
if OPT_USEETF
    isMicro = isMicro(2:end,:);
end

% End-of-Month
signals_LF = signals_LF(1:end-OPT_LAG,:,:);
signals_HF = signals_HF(1:end-OPT_LAG,:,:);
isMicro    = isMicro(1:end-OPT_LAG,:);
% cap        = cap(1:end-OPT_LAG,:);
inan_sig   = inan_sig(1:end-OPT_LAG,:,:);

% Lag forward
hpr   = hpr(1+OPT_LAG:end,:);
rf    = rf(1+OPT_LAG:end,:);
mdate = mdate(1+OPT_LAG:end,:);
idec  = mdate > 200104;

% Filter micro
hpr(isMicro) = NaN;
%% BAB ret
[ptfret, avgsig,wl,wh] = deal(cell(nsig,2));
for ii = 1:nsig
    [ptfret{ii,1},~,wl{ii,1},wh{ii,1},avgsig{ii,1}] = bab(hpr,signals_LF(:,:,ii),rf);
    [ptfret{ii,2},~,wl{ii,2},wh{ii,2},avgsig{ii,2}] = bab(hpr,signals_HF(:,:,ii),rf);
end
%% Desc stats

% Correlation
snames           = {'babm','babq','babs','baby','rbabm','rbabq','rbabs','rbaby'};
order            = reshape(reshape(1:nsig*2,nsig,2)',1,nsig*2);
allsig           = cat(3,signals_LF,signals_HF);
correlations     = corrxs(allsig(:,:,order),snames(order));
correlations_dec = corrxs(allsig(idec,:,order),snames(order));

% Plot
figure
dt = serial2datetime(datenum(1993,(1:size(signals_LF,1))+2,1)-1);
for ii = 1:nsig*2
    subplot(nsig,2,ii)
    inan      = isnan(ptfret{order(ii)});
    lvl       = cumprod(1+nan2zero(ptfret{order(ii)}));
    lvl(inan) = NaN;
    plot(dt,lvl)
    set(gca, 'TickLabelInterpreter','latex','Layer','Top')
end

% Plot
figure
dt = serial2datetime(datenum(1993,(1:size(signals_LF,1))+2,1)-1);
X  = datenum([min(dt(idec)), min(dt(idec)), max(dt(idec)),max(dt(idec))]);
for ii = 1:nsig*2
    subplot(nsig,2,ii)
    inan      = isnan(ptfret{order(ii)}(idec,:));
    lvl       = cumprod(1+nan2zero(ptfret{order(ii)}(idec,:)));
    lvl(inan) = NaN;
    plot(dt(idec),lvl)
    set(gca, 'TickLabelInterpreter','latex','Layer','Top')
end

% Desc stats
desc   = cellfun(@(x) stratstats(dt,x*100,'Frequency','m','IsPercentageReturn',true), ptfret,'un',0);
desc   = cellfun(@(x) renameVarNames(x',{'Low','High','BAB'}), desc,'un',0);
catfun = @(sig,stats) [stats; array2table(nanmean(sig),'VariableNames',{'Low','High','BAB'},'RowNames',{'Avgsig'})];
desc   = cellfun(catfun,avgsig,desc,'un',0);
out    = cell2mat(cellfun(@(x) x{:,:}, desc,'un',0));

% Desc stats, post decimalization
% Jan 2001 decimalization in NYSE and April in NASDAQ, take from May
desc2  = cellfun(@(x) stratstats(dt(idec),x(idec,:)*100,'Frequency','m','IsPercentageReturn',true), ptfret,'un',0);
desc2  = cellfun(@(x) renameVarNames(x',{'Low','High','BAB'}), desc2,'un',0);
catfun = @(sig,stats) [stats; array2table(nanmean(sig(idec,:)),'VariableNames',{'Low','High','BAB'},'RowNames',{'Avgsig'})];
desc2  = cellfun(catfun,avgsig,desc2,'un',0);
out    = cell2mat(cellfun(@(x) x{:,:}, desc2,'un',0));
%% Tests
[~,pValST,Zjk]   = cellfun(@(high,low,M) sharpetest(high(:,end)   , low(:,end)   ,M,[],'both'), ptfret(:,2),ptfret(:,1), {1,3,6,12}');
[~,pValST2,Zjk2] = cellfun(@(high,low,M) sharpetest(high(idec,end), low(idec,end),M,[],'both'), ptfret(:,2),ptfret(:,1), {1,3,6,12}');
% 
% [coeff,se,tratio,pval]     = cellfun(@(high,low) regressHighOnLow(high(:,end), low(:,end)), ptfret(:,2),ptfret(:,1),'un',0);
% [coeff2,se2,tratio2,pval2] = cellfun(@(high,low) regressHighOnLow(high(idec,end), low(idec,end)), ptfret(:,2),ptfret(:,1),'un',0);

% [se3, pval3] = cellfun(@(high,low) sharpeHAC([high(:,end), low(:,end)]), ptfret(:,2),ptfret(:,1));
% [se4, pval4] = cellfun(@(high,low) sharpeHAC([high(idec,end), low(idec,end)]), ptfret(:,2),ptfret(:,1));

% Ledoit, Wolf (2008) boot-strapped test
rng default
% OPT_BLOCKS_DEC = cellfun(@(high,low) blockSizeCalibrate([high(:,end), low(:,end)]), ptfret(:,2),ptfret(:,1));
% OPT_BLOCKS_DEC = cellfun(@(high,low) blockSizeCalibrate([high(idec,end), low(idec,end)]), ptfret(:,2),ptfret(:,1));
rng default
pval5 = cellfun(@(high,low,b) bootInference([high(idec,end), low(idec,end)],b,[],[],0), ptfret(:,2),ptfret(:,1),{2 2 2 2}');
rng default
pval6 = cellfun(@(high,low,b) bootInference([high(~isnan(high(:,end)),end), low(~isnan(low(:,end)),end)],b,[],[],0), ptfret(:,2),ptfret(:,1),{2 2 2 2}');

%% Conditioning: illiquidity
illiq = loadresults('illiq');
illiq = illiq(1:end-OPT_LAG,:,:);

SR = repmat({NaN(OPT_PTF_DB,2)},nsig,1);

for ii = 1:nsig
    il                   = illiq;
    il(inan_sig(:,:,ii)) = NaN;
    [bins,counts]        = binSignal(il,'PortfolioNumber',OPT_PTF_DB);

    for jj = 1:OPT_PTF_DB
        iLiq = bins == jj;

        % Low freq
        signal                = signals_LF(:,:,ii);
        signal(~iLiq)         = NaN;
        [tmpret,~,~,~,tmpsig] = bab(hpr,signal,rf);

        s = stratstats(dt(idec), tmpret(idec,:)*100,'Frequency','m','IsPercentageReturn',true);

        SR{ii}(jj,1) = s.SR(end);

        % High freq
        signal                = signals_HF(:,:,ii);
        signal(~iLiq)         = NaN;
        [tmpret,~,~,~,tmpsig] = bab(hpr,signal,rf);

        s            = stratstats(dt(idec), tmpret(idec,:)*100,'Frequency','m','IsPercentageReturn',true);
        SR{ii}(jj,2) = s.SR(end);
    end
end
celldisp(SR)
%% Conditioning: mkt cap
cap   = getMktCap(xstr2num(permno),[],1);
idx   = ismember(cap.Date, date);
cap   = cap(idx,:);

[~,pos] = unique(cap.Date/100,'last');
cap     = cap{pos,2:end};
cap     = cap(1:end-OPT_LAG,:,:);

% Check correlation with signals
% Note: illiquidity is -90% corr with log(size)
corrxs(cat(3,allsig(idec,:,order),illiq(idec,:), log(cap(idec,:))),[snames(order), 'illiq','cap']);

SR = repmat({NaN(OPT_PTF_DB,2)},nsig,1);

for ii = 1:nsig
    c                   = cap;
    c(inan_sig(:,:,ii)) = NaN;
    [bins,counts]       = binSignal(c,'PortfolioNumber',OPT_PTF_DB);

    for jj = 1:OPT_PTF_DB
        iCap = bins == jj;

        % Low freq
        signal                = signals_LF(:,:,ii);
        signal(~iCap)         = NaN;
        [tmpret,~,~,~,tmpsig] = bab(hpr,signal,rf);

        s = stratstats(dt(idec), tmpret(idec,:)*100,'Frequency','m','IsPercentageReturn',true);

        SR{ii}(jj,1) = s.SR(end);

        % High freq
        signal                = signals_HF(:,:,ii);
        signal(~iCap)         = NaN;
        [tmpret,~,~,~,tmpsig] = bab(hpr,signal,rf);

        s            = stratstats(dt(idec), tmpret(idec,:)*100,'Frequency','m','IsPercentageReturn',true);
        SR{ii}(jj,2) = s.SR(end);
    end
end
celldisp(SR)
%% Risk-adjustment
factors = loadresults('RAfactors');
factors = factors(ismember(factors.Date, mdate),:);
%% Profitability
[idx,msubs] = ismember(date/100,mdate);
msubs = msubs(idx);

% Apply filters inherited from returns and exclude hpr NaNs
inan = isnan(hpr);
hprd = ret(idx,:);
hprd(inan(msubs,:)) = NaN;

% Loop for all signals
% [num,n] = deal(zeros(opt.Ptf_num_univ,20,nsignals));
for jj = 1:nsig
    rh     = nansum(hprd .* wh{jj,1}(msubs,:),2);
    rl     = nansum(hprd .* wl{jj,1}(msubs,:),2);
    bh     = nansum(signals_LF(:,:,jj) .* wh{jj,1},2);
    bl     = nansum(signals_LF(:,:,jj) .* wl{jj,1},2);
    ptfret = (rl-ff.RF(idx)/100)./bl(msubs,:) - (rh-ff.RF(idx)/100)./bh(msubs,:);
    
	rh     = nansum(hprd .* wh{jj,2}(msubs,:),2);
    rl     = nansum(hprd .* wl{jj,2}(msubs,:),2);
    bh      = nansum(signals_HF(:,:,jj) .* wh{jj,2},2);
    bl      = nansum(signals_HF(:,:,jj) .* wl{jj,2},2);
    ptfret  = [ptfret, (rl-ff.RF(idx)/100)./bl(msubs,:) - (rh-ff.RF(idx)/100)./bh(msubs,:)];
    
    ptfret(isinf(ptfret)) = NaN;
    figure
    plot(cumprod(nan2zero(ptfret)+1))
%     
%     % Series to keep
%     ikeep = ~inan & bin ~= 0;
%     
%     for ii = 1:max(msubs)
%         imonth   = msubs == ii;
%         retslice = ret(imonth, ikeep(ii,:));
%         if ~isempty(retslice)
%             % Accumulation indexes
%             [row,col] = ndgrid(1:size(retslice,1), bin(ii, ikeep(ii,:)));
%             subs      = [row(:), col(:)];
%         
%             % Mean as S(sum)/S(n)
%             tmp         = accumarray(subs, retslice(:));
%             num(:,:,jj) = num(:,:,jj) + tmp([1:10,end-9:end],:)';
%             tmp         = accumarray(subs, 1);
%             n(:,:,jj)   = n(:,:,jj) + tmp([1:10,end-9:end],:)';
%         end
%     end
end
