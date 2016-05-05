%% Options
OPT_FREQ = 75;

OPT_LAG    = 1;
OPT_PTF_UN = 5;
OPT_PTF_DB = 5;

OPT_SHRINK = [0.4,0.6,0.6,0.6];
% OPT_SHRINK = [1,1,1,1];

% OPT_BLOCKS_DEC = {4 2 1 1}'; % whole horizon
OPT_BLOCKS_DEC = {1 6 4 2}';
% OPT_BLOCKS_DEC = {6;2;8;8} % Preav
%% Data
load(sprintf('results\\alldata_beta%d',OPT_FREQ))

%% Signals
[signals_LF, hpr, rf, mdate] = make_signals_LF(ret,date,ff,OPT_SHRINK);
signals_HF                   = make_signals_HF(xstr2num(permno),date,beta,OPT_SHRINK);

% signals_LF = filter_signal(arima(1,0,1), signals_LF);
% signals_HF = filter_signal(arima(1,0,1), signals_HF);

% NaN-out intersection
nsig     = size(signals_LF,3);
inan_sig = false(size(signals_LF));
for ii = 1:nsig
    inan_sig(:,:,ii) = isnan(signals_LF(:,:,ii)) | isnan(signals_HF(:,:,ii));
end
signals_LF(inan_sig) = NaN;
signals_HF(inan_sig) = NaN;

%% Lag

% End-of-Month
signals_LF = signals_LF(1:end-OPT_LAG,:,:);
signals_HF = signals_HF(1:end-OPT_LAG,:,:);
isMicro    = isMicro(1:end-OPT_LAG,:);
inan_sig   = inan_sig(1:end-OPT_LAG,:,:);

% Lag forward
hpr   = hpr(1+OPT_LAG:end,:);
rf    = rf(1+OPT_LAG:end,:);
mdate = mdate(1+OPT_LAG:end,:);
idec  = mdate > 200104;
dt    = yyyymmdd2datetime(mdate*100);

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
outorder         = reshape(reshape(1:nsig*2,nsig,2)',1,nsig*2);
allsig           = cat(3,signals_LF,signals_HF);
correlations     = corrxs(allsig(:,:,outorder),snames(outorder));
correlations_dec = corrxs(allsig(idec,:,outorder),snames(outorder));

% % Plot
% figure
% for ii = 1:nsig*2
%     subplot(nsig,2,ii)
%     inan      = isnan(ptfret{outorder(ii)});
%     lvl       = cumprod(1+nan2zero(ptfret{outorder(ii)}));
%     lvl(inan) = NaN;
%     plot(dt,lvl)
%     set(gca, 'TickLabelInterpreter','latex','Layer','Top')
% end
% 
% % Plot
% figure
% for ii = 1:nsig*2
%     subplot(nsig,2,ii)
%     inan      = isnan(ptfret{outorder(ii)}(idec,:));
%     lvl       = cumprod(1+nan2zero(ptfret{outorder(ii)}(idec,:)));
%     lvl(inan) = NaN;
%     plot(dt(idec),lvl)
%     set(gca, 'TickLabelInterpreter','latex','Layer','Top')
% end

% Desc stats
desc   = cellfun(@(x) stratstats(dt,x*100,'Frequency','m','IsPercentageReturn',true), ptfret,'un',0);
desc   = cellfun(@(x) renameVarNames(x',{'Low','High','BAB'}), desc,'un',0);
catfun = @(sig,stats) [stats; array2table(nanmean(sig),'VariableNames',{'Low','High','BAB'},'RowNames',{'Avgsig'})];
desc   = array2table(cellfun(catfun,avgsig,desc,'un',0),'VariableNames',{'LF','HF'},'RowNames',{'m','s','q','y'});

% Desc stats, post decimalization
% Jan 2001 decimalization in NYSE and April in NASDAQ, take from May
desc2  = cellfun(@(x) stratstats(dt(idec),x(idec,:)*100,'Frequency','m','IsPercentageReturn',true), ptfret,'un',0);
desc2  = cellfun(@(x) renameVarNames(x',{'Low','High','BAB'}), desc2,'un',0);
catfun = @(sig,stats) [stats; array2table(nanmean(sig(idec,:)),'VariableNames',{'Low','High','BAB'},'RowNames',{'Avgsig'})];
desc2  = array2table(cellfun(catfun,avgsig,desc2,'un',0),'VariableNames',{'LF','HF'},'RowNames',{'m','s','q','y'});

SR = NaN(4,2);
for ii = 1:4
    SR(ii,:) = [desc2.LF{ii}{'SR','BAB'} desc2.HF{ii}{'SR','BAB'}];
end

%% Tests
% [~,pValST,Zjk]   = cellfun(@(high,low,M) sharpetest(high(:,end)   , low(:,end)   ,M,[],'both'), ptfret(:,2),ptfret(:,1), {1,3,6,12}');
[~,pValST2,Zjk2] = cellfun(@(high,low,M) sharpetest(high(idec,end), low(idec,end),M,[],'both'), ptfret(:,2),ptfret(:,1), {1,3,6,12}');
% 
% [coeff,se,tratio,pval]     = cellfun(@(high,low) regressHighOnLow(high(:,end), low(:,end)), ptfret(:,2),ptfret(:,1),'un',0);
% [coeff2,se2,tratio2,pval2] = cellfun(@(high,low) regressHighOnLow(high(idec,end), low(idec,end)), ptfret(:,2),ptfret(:,1),'un',0);

% [se3, pval3] = cellfun(@(high,low) sharpeHAC([high(:,end), low(:,end)]), ptfret(:,2),ptfret(:,1));
% [se4, pval4] = cellfun(@(high,low) sharpeHAC([high(idec,end), low(idec,end)]), ptfret(:,2),ptfret(:,1));

% Ledoit, Wolf (2008) boot-strapped test
% rng default
% OPT_BLOCKS_DEC = cellfun(@(high,low) blockSizeCalibrate([high(:,end), low(:,end)]), ptfret(:,2),ptfret(:,1));
% OPT_BLOCKS_DEC = cellfun(@(high,low) blockSizeCalibrate([high(idec,end), low(idec,end)]), ptfret(:,2),ptfret(:,1));
% rng default
% pval5 = cellfun(@(high,low,b) bootInference([high(~isnan(high(:,end)),end), low(~isnan(low(:,end)),end)],b,[],[],0), ptfret(:,2),ptfret(:,1),OPT_BLOCKS_DEC);
rng default
pval6 = cellfun(@(high,low,b) bootInference([high(idec,end), low(idec,end)],b,[],[],0), ptfret(:,2),ptfret(:,1),OPT_BLOCKS_DEC);
%% Conditioning: illiquidity
illiq = loadresults('illiq');
illiq = illiq(1:end-OPT_LAG,:,:);

M  = {1,3,6,12};
SR = repmat({NaN(OPT_PTF_DB,2)},nsig,1);
[pvalill1, pvalill2] = deal(repmat({NaN(OPT_PTF_DB,1)},nsig,1));
for ii = 1:nsig
    il                   = illiq;
    il(inan_sig(:,:,ii)) = NaN;
    [bins,counts]        = binSignal(il,'PortfolioNumber',OPT_PTF_DB);

    for jj = 1:OPT_PTF_DB
        iLiq = bins == jj;

        % Low freq
        signal                = signals_LF(:,:,ii);
        signal(~iLiq)         = NaN;
        [retlf,~,~,~,tmpsig] = bab(hpr,signal,rf);

        s = stratstats(dt(idec), retlf(idec,:)*100,'Frequency','m','IsPercentageReturn',true);
        SR{ii}(jj,1) = s.SR(end);

        % High freq
        signal                = signals_HF(:,:,ii);
        signal(~iLiq)         = NaN;
        [rethf,~,~,~,tmpsig] = bab(hpr,signal,rf);

        s            = stratstats(dt(idec), rethf(idec,:)*100,'Frequency','m','IsPercentageReturn',true);
        SR{ii}(jj,2) = s.SR(end);
        
        [~,pvalill1{ii}(jj)] = sharpetest(rethf(idec,end), retlf(idec,end),M{ii},[],'both');
%         pvalill2{ii}(jj) = bootInference([rethf(idec,end), retlf(idec,end)],OPT_BLOCKS_DEC{ii},[],[],0);
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
corrxs(cat(3,allsig(idec,:,outorder),illiq(idec,:), log(cap(idec,:))),[snames(outorder), 'illiq','cap']);

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
    ptfretd = (rl-ff.RF(idx)/100)./bl(msubs,:) - (rh-ff.RF(idx)/100)./bh(msubs,:);
    
	rh     = nansum(hprd .* wh{jj,2}(msubs,:),2);
    rl     = nansum(hprd .* wl{jj,2}(msubs,:),2);
    bh      = nansum(signals_HF(:,:,jj) .* wh{jj,2},2);
    bl      = nansum(signals_HF(:,:,jj) .* wl{jj,2},2);
    ptfretd  = [ptfretd, (rl-ff.RF(idx)/100)./bl(msubs,:) - (rh-ff.RF(idx)/100)./bh(msubs,:)];
    
    ptfretd(isinf(ptfretd)) = NaN;
    figure
    plot(cumprod(nan2zero(ptfretd)+1))

    retbyday(date(idx), ptfretd)
end
idt = date(idx) > 20010431;
stratstats(date(date>20010431),ptfretd(idt,:)*100,'Frequency','d','IsPercentageReturn',true)
