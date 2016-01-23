function [stratret, stratlvl] = zeroptf(tb, score)
% [stratret, tbstats, tbarets] = zeroptf(tb)
%
%   Where TB should have (in no order):  
%       ID | Date | Score | Ret

warning off MATLAB:table:ModifiedVarnames
% Unstack returns
ret    = unstack(tb(:,{'ID','Date','Ret'}),'Ret','ID');
ret    = sortrows(ret,'Date');
dates  = uint32(ret.Date);

% Unstack scores (eventually)
if nargin < 2
    score = unstack(tb(:,{'ID','Date','Score'}),'Score','ID');
    score = sortrows(score,'Date');
elseif ~isequal(score.Properties.VariableNames,ret.Properties.VariableNames) ||...
       ~isequal(score.Date, ret.Date)
    error('zeroptf:wrongScores','Check that SCORE has same IDs and dates.')
end
warning on MATLAB:table:ModifiedVarnames

ret   = table2array(ret(:,2:end));
score = table2array(score(:,2:end));

% Remove rows with no score
inan  = all(isnan(score),2);
score = score(~inan,:);
ret   = ret(~inan,:);
dates = dates(~inan);

% Beginning of month rebalancing scheme
[~, ~, subs] = unique(dates./100);
rebdate      = find([false; diff(subs)>0]);

N        = numel(dates);
stratret = table(dates, NaN(N,1),NaN(N,1),NaN(N,1),...
                 'VariableNames',{'Date','Top','Bottom','Ptf'});
for ii = 1:numel(rebdate)
    [iBottom, iTop] = deal(false(1,size(ret,2)));
    % Alive on rebalancing date
    ifut    = subs == ii+1;
    isalive = ~isnan(score(rebdate(ii),:));
    
    % Score ranking
    scores           = score(rebdate(ii)-1,isalive);
    ptiles           = prctile(scores,[10,90]);
    iBottom(isalive) = scores <= ptiles(1);
    iTop (isalive)   = scores >= ptiles(2);
       
    % Equal weighted strategy
    stratret.Top(ifut)    = nanmean(ret(ifut,iTop),2);
    stratret.Bottom(ifut) = nanmean(ret(ifut,iBottom),2);
end
stratret.Ptf = stratret.Top-stratret.Bottom;

if nargout == 2
    from            = find(~isnan(stratret.Ptf),1,'first');
    fun             = @(x) [NaN(from-2,1); cumprod([1; x(from:end)+1])];
    stratlvl        = stratret;
    stratlvl(:,2:4) = varfun(fun, stratlvl, 'InputVariables',2:4);
end
end