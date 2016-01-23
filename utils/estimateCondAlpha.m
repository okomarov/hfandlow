function scores = estimateCondAlpha(lookback, rebdates, betas, rets, f, C)
% estimateCondAlpha(betas, rets)
% Sorts betas by Date and Permno

% Cond alphas
% From r_it = alpha_{i,t-1} + f_t * beta_{i,t-1} + e_it:
% 1) estimate the daily risk-adjusted return ^Z_it = r_it - f_t * ^beta_{i,t-1}
% 2) estimate the conditional alpha E[^Z_it | C_{t-1}]

% Unstack returns
warning off MATLAB:table:ModifiedVarnames
rets = rets(~isnan(rets.RetCC),{'Date','Permno','RetCC'});
rets = unstack(rets,'RetCC','Permno');
rets = sortrows(rets,'Date');
warning on MATLAB:table:ModifiedVarnames

% Intersect columns
[~,ia,ib] = intersect(getVariableNames(betas),getVariableNames(rets));
betas     = betas(:,ia);
rets      = rets(:,ib);

% Ensure f sorted by date
f = sortrows(f(:,{'Date','RetCC'}),'Date');

% Intersect dates
refdates = intersect(rets.Date, f.Date);

% Time align: r_t = f_t * BETA_{t-1}' and E[^Z_it | C_{t-1}]
notrail = true;
% Sample at t-1
tminus1 = serial2yyyymmdd(yyyymmdd2serial(refdates) - 1);
betas   = sampledates(betas, tminus1,  notrail);
C       = sampledates(C    , tminus1,  notrail);
% Sample at t
rets    = sampledates(rets,  refdates, notrail);
f       = sampledates(f   ,  refdates, notrail);

% Position of scoring days, rebalance date - 1
[~,rebpos] = ismember(rebdates,refdates);
scorepos   = rebpos - 1;
scorepos   = scorepos(scorepos >= lookback);

% Extract data
vnames = getVariableNames(rets);
Date   = rets.Date;
C      = C{:,2:end};
rets   = rets{:,2:end};
betas  = betas{:,2:end};
f      = f{:,2:end};


[nobs,nseries] = size(rets);
scores         = NaN(nobs,nseries);
for row = scorepos(:)' 
    pos = row-lookback+1:row;
    % Residual returns ^Z_it = r_it - f_t * ^beta_{i,t-1}
    Z   = rets(pos,:) - bsxfun(@times, f(pos), betas(pos,:));
    for col = 1:nseries
        y = Z(:,col);
        if nnz(~isnan(y)) > lookback * 0.5
            % conditional alpha E[^Z_it | C_{t-1}]
            mdl             = regstats(y,C(pos,:),'linear','beta');
            scores(row,col) = mdl.beta(1);
        end
    end
end
scores = nanfillts(scores);
scores = [table(Date), array2table(scores,'VariableNames',vnames(2:end))];
end
