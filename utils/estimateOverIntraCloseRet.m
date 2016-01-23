function res = estimateOverIntraCloseRet
% Calulate overnight return from CRSP's close-to-close return and
% TAQ's open-to-close (intraday) return.


% Load Permno - Date pairs
mst = loadresults('masterPermno','..\results');
mst = addPermno(mst);

% Open-to-close returns
fl        = loadresults('sampleFirstLast','..\results');
retOC     = fl.LastPrice./fl.FirstPrice-1;
[~,pos]   = ismembIdDate(mst.Id, mst.Date, fl.Id, fl.Date);
mst.RetOC = retOC(pos);

% Load dsfquery
dsf = loadresults('dsfquery','..\results');

% % Winsorize returns at 0.1 and 99.9%
% ptiles = prctile(dsfquery.Ret,[0.1,99.9]);
% % boxplot(dsfquery.Ret)
% % idx    = in(dsfquery.Ret, ptiles);
% % boxplot(dsfquery.Ret(idx))
% dsfquery.Ret(dsfquery.Ret < ptiles(1)) = ptiles(1);
% dsfquery.Ret(dsfquery.Ret > ptiles(2)) = ptiles(2);

% Close-to-close return
[idx,pos]      = ismembIdDate(mst.Permno, mst.Date, dsf.Permno, dsf.Date);
mst.RetCC      = NaN(size(mst.Permno));
mst.RetCC(idx) = dsf.Ret(pos(idx));

% Overnight return
mst.RetCO = log((1 + mst.RetCC)./(1 + mst.RetOC ));

res      = mst;
filename = sprintf('%s_%s.mat',datestr(now,'yyyymmdd_HHMM'),'return_intraday_overnight');
save(fullfile('.\results\',filename), 'res')
end