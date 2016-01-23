function betas = estimateBetaComponents(freq, useon, useproxy, grid)
if nargin < 1 || isempty(freq),     freq     = 5;     end
if nargin < 2 || isempty(useon),    useon    = true;  end
if nargin < 3 || isempty(useproxy), useproxy = false; end
if nargin < 4 || isempty(grid),     grid     = [];    end


writeto = '.\results\';

try
    fprintf('%s: loading betacomponents at %d min.\n', mfilename, freq)
    name  = matname('betacomponents',freq, useon, useproxy);
    betas = loadresults(name);
catch
    fprintf('%s: betacomponents not found. Estimating...\n', mfilename)

    % If sampling freq is a multiple of 5, use 5min sampled data and
    % compose later
    isComposable = mod(freq,5) == 0 && nargin == 4;
    if isComposable
        path2data = '..\data\TAQ\sampled\5min\nobad_vw';        
    else
        path2data = fullfile('..\data\TAQ\sampled\', sprintf('%dmin', freq),'nobad_vw');
    end

    % Eventully Sample data
    if ~isComposable && (exist(path2data,'dir') ~= 7 || numel(dir(path2data)) <= 2)
        fprintf('%s: sampling data at %d min.\n', mfilename, freq)
        step    = freq/(60*24);
        grid    = (9.5/24:step:16/24)';
        fmtname = sprintf('S%dm_%%04d.mat',freq);
        sampleData(grid, path2data, fmtname);
    end
    
    if isComposable
        [spret, reton] = getMktRet(5, useproxy, useon, path2data);
    else
        [spret, reton] = getMktRet(freq, useproxy, useon, path2data);
    end

    load(fullfile(path2data,'master'),'-mat','mst');
    cached = cacheReplicateSpret(spret,mst);

    % Cache overnight returns
    % Note: the match is on Date -Id rather than Date - Permno to avoid the
    %       duplication of overnight returns that comes from the
    %       Date - Permno mapping to Date - Permno/Id
    if useon
        fprintf('%s: adding overnight returns to all other series.\n', mfilename)
        reton.File      = zeros(size(reton,1),1,'uint16');
        [idx,pos]       = ismembIdDate(reton.Id,reton.Date, mst.Id, mst.Date);
        reton.File(idx) = mst.File(pos(idx));
        reton           = reton(idx,{'Permno','Date','RetCO','File'});
        cached          = [cached,...
                           accumarray(reton.File,(1:size(reton))',[],@(x) {reton(x,{'Permno','Date','RetCO'})})];
    end
    clear mst reton pos idx

    % Calculate beta components: sum(r*benchr) and sum(benchr^2)
    fprintf('%s: creating betacomponents at %d min.\n', mfilename, freq)
    [betas,filename] = AnalyzeHflow('betacomponents', [], cached, path2data,[],8,grid);

    % Rename to append the sampling frequency
    name        = regexp(filename,'\w+?(?=\.mat)','match','once');
    name        = [matname(name,freq, useon, useproxy),'.mat'];
    newfullname = fullfile(writeto, name);
    movefile(fullfile(writeto,filename), newfullname);
end
end

function [mktret,retOn] = getMktRet(freq,useproxy,useon,path2data)
% Make intraday sp500
if useproxy
    name = sprintf('sp500proxy%dm',freq);
    try
        mktret = loadresults(name);
    catch
        fprintf('%s: creating ssp500proxy at %d min.\n', mfilename, freq)
        mktret = sp500intraday(path2data);
    end

    % ETF Spyders
else
    mktret = getSpy(freq);
end

mktret = [mktret.Datetime, [NaN; diff(log(mktret.Price))]];

% Deal with overnight
iNotOn = [false; diff(rem(mktret(:,1),1)) >= 0];
if useon
    fprintf('%s: adding overnight returns to the index.\n', mfilename)
    posOn                = find(~iNotOn);
    retOn                = loadresults('return_intraday_overnight');
    spreton              = retOn(retOn.Permno == 84398,:);
    [idx,pos]            = ismember(serial2yyyymmdd(mktret(posOn,1)),spreton.Date);
    mktret(posOn(idx),2) = spreton.RetCO(pos(idx));
else
    retOn  = [];
    mktret = mktret(iNotOn,:);
end
end

function cached = cacheReplicateSpret(spret, mst)
% cached = {cellarray of returns, cellarray of dates} - [Nfiles by 2]

fprintf('%s: caching index returns by days.\n', mfilename)
nfiles = max(mst.File);
cached = cell(nfiles,2);

% Cache SP500 returns by days
dates           = fix(spret(:,1));
[spdays,~,subs] = unique(dates,'stable');
spret           = cache2cell(spret, subs);

% Replicate to match master
unMstDates = accumarray(mst.File, mst.Date,[],@(x){yyyymmdd2serial(unique(x))});
for ii = 1:nfiles
    [~,pos]      = ismember(unMstDates{ii}, spdays);
    nnzero       = pos ~= 0;
    isp          = ismember(spdays, unMstDates{ii});
    cached(ii,:) = {spret(isp,:) spdays(pos(nnzero))};
end
end

function name = matname(name, freq, useon, useproxy)
if useon,    useon    = 'on';   else useon    = ''; end
if useproxy, useproxy = 'prx';  else useproxy = ''; end
name = sprintf('%s%dm%s%s', name, freq, useon, useproxy);
end