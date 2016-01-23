function [coeff, se, tratio, pval] = risk_adjust(retcell, factors, ptfnum)

% Preallocate
N          = numel(retcell);
[coeff,se] = deal(cell(N));

% (Long - short)_t = a + X_t * b
for ii = 1:N
    r = retcell{ii};

    % Percentage long - short ptf return
    ptfret             = (r(:,end) - r(:,1))*100;
    [coeff{ii},se{ii}] = regressions(ptfret,factors);
end
coeff = cat(1,coeff{:});
se    = cat(1,se{:});

% tstat and pvalues
tratio = coeff./se;
pval   = 2 * normcdf(-abs(tratio));
end

function [coeff,se] = regressions(y, ff)

% Calculates regression coefficients and se, and stores them in the
% following format:
%
%     Model 1 | Model 2 | ...
% c
% x1
% x2
% ...
%
% R^2

% filter out nans
nonan = ~isnan(y);
y     = y(nonan);
ff    = ff(nonan,:);
nobs  = nnz(nonan);
l     = ones(nobs,1);

% Preallocations
ntests     = 6;
nfactors   = 13 + 1;
[coeff,se] = deal(NaN(nfactors,ntests));

opts = {'intercept',false,'display','off','type','HAC','bandwidth',floor(4*(nobs/100)^(2/9))+1,'weights','BT'};
f    = @(x)     hac(x, y, opts{:});
rsq  = @(X,b)   var(X * b)/var(y);
adjr = @(X,b,p) 1 - (1-rsq(X,b))*(nobs-1)/(nobs-p-1);

% Excess
col                            = 1;
row                            = 1;
[~,se(row,col),coeff(row,col)] = f(l); col = col + 1;
% FF3
row                            = 1:4;
X                              = [l ff{:,{'MktMinusRF', 'SMB', 'HML'}}];
[~,se(row,col),coeff(row,col)] = f(X);
coeff(end-1,col)               = adjr(X,coeff(row,col), size(X,2)-1); col = col + 1;
% HXZ'15
row                            = [1:3,7:8];
X                              = [l ff{:,{'MKT', 'ME', 'IA','ROE'}}];
[~,se(row,col),coeff(row,col)] = f(X);
coeff(end-1,col)               = adjr(X,coeff(row,col), size(X,2)-1); col = col + 1;
% FF5
row                            = 1:6;
X                              = [l ff{:,{'MktMinusRF', 'SMB', 'HML', 'RMW', 'CMA'}}];
[~,se(row,col),coeff(row,col)] = f(X);
coeff(end-1,col)               = adjr(X,coeff(row,col), size(X,2)-1); col = col + 1;
% FF5 + LIQ (Pastor)
row                            = [1:6,9];
X                              = [l ff{:,{'MktMinusRF', 'SMB', 'HML', 'RMW', 'CMA','LIQ_p'}}];
[~,se(row,col),coeff(row,col)] = f(X);
coeff(end-1,col)               = adjr(X,coeff(row,col), size(X,2)-1); col = col + 1;
% FF5 + MOM + TSMOM + STREV
row                            = [1:6,10:12];
X                              = [l ff{:,{'MktMinusRF', 'SMB', 'HML', 'RMW', 'CMA','Mom','TSmom','ST_Rev'}}];
[~,se(row,col),coeff(row,col)] = f(X);
coeff(end-1,col)               = adjr(X,coeff(row,col), size(X,2)-1);
end