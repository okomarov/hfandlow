function betas = estimateBetaComponents_refresh(useon, badPriceMult, consolidateType, minRefStep, method)
if nargin < 1 || isempty(useon), useon           = true;                end
if nargin < 2,                   badPriceMult    = [];                  end
if nargin < 3,                   consolidateType = 'volumeWeighted';    end
if nargin < 4,                   minRefStep      = 30/(24*60);          end
if nargin < 5,                   method          = 'refresh';           end

try
    fprintf('%s: loading betacomponents_%s.\n', mfilename,method)
    betas = loadresults(['betacomponents_', method]);
catch
    fprintf('%s: betacomponents not found. Estimating...\n', mfilename)

    path2data = '..\data\TAQ\';

    mst = prepareMst();

    % Cache overnight returns
    % Note: the match is on Date -Id rather than Date - Permno to avoid the
    %       duplication of overnight returns that comes from the
    %       Date - Permno mapping to Date - Permno/Id
    if useon
        fprintf('%s: adding overnight returns to all other series.\n', mfilename)
        reton            = loadresults('return_intraday_overnight');
        reton.File       = zeros(size(reton,1),1,'uint16');
        [idx,pos]        = ismembIdDate(mst.Id,mst.Date, reton.Id, reton.Date);
        mst.RetCO(idx,1) = reton.RetCO(pos(idx));
    end
    clear reton pos idx

    spymst                                     = mst(mst.Permno == 84398,:);
    spymst(:,{'Permno','MedPrice','Isbadday'}) = [];
    mst(:,{'From','To'})                       = [];
    mst                                        = cache2cell(mst,mst.File);
    spymst                                     = cacheReplicateSpret(spymst, mst);

    % Calculate beta components: sum(r*benchr) and sum(benchr^2)
    fprintf('%s: creating betacomponents_%s.\n', mfilename,method)
    opt   = struct('HasOvernight',useon,'BadPriceMultiplier',badPriceMult,...
                   'TimestampConsolidation',consolidateType,'MinRefreshStep',minRefStep);

    switch method
        case 'preaverage'
            betas = Analyze('betacomponents_preav', [], [mst,spymst], path2data,[],6,opt);
        case 'refresh'
            betas = Analyze('betacomponents_refresh', [], [mst,spymst], path2data,[],6,opt);
    end
end
end

function cached = cacheReplicateSpret(spymst, mst)
% cached = {cellarray of returns, cellarray of dates} - [Nfiles by 2]

fprintf('%s: caching index returns by days.\n', mfilename)
nfiles = numel(mst);
cached = cell(nfiles,2);

% Replicate to match master
for ii = 1:nfiles
    unMstDates   = unique(mst{ii}.Date);
    [idx,pos]    = ismember(unMstDates, spymst.Date);
    cached(ii,:) = {unMstDates(idx), spymst(pos(idx),:)};
end
end
