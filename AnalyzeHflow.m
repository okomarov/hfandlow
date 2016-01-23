function [res, filename] = AnalyzeHflow(fun, varnames, cached, path2data, debug, poolcores, varargin)
% ANALYZE Executes specified fun in parallel on the whole database (all .mat files)
%
%   ANALYZE(FUN, VARNAMES) FUN should a string with the name of one of
%                          the following sub-functions:
%                               - 'dailystats'
%                               - 'badprices'
%                               - 'avgtimestep'
%                          VARNAMES is a cell-array of strings (or string)
%                          with the VarNames of the dataset with the results
%
%   ANALYZE(..., PATH2DATA) If you wanna use other than '.\data\TAQ\T*.mat'
%                           files (default), then specify a different
%                           PATH2DATA with the initial pattern of the name,
%                           e.g '.\data\TAQ\sampled\S5m_*.mat'
%
%   ANALYZE(..., CACHED) Some FUN might require pre-cached results which where
%                        run on the whole database.
%                        Check the specific sub-function for the format of the
%                        needed CACHED results.
%   ANALYZE(..., DEBUG) Run execution sequentially, i.e. not in parallel, to be
%                       able to step through the code in debug mode.
if nargin < 2 || isempty(varnames);  varnames  = {'data','mst','ids'};  end
if nargin < 3,                       cached    = [];                    end
if nargin < 4 || isempty(path2data); path2data = '.\data\TAQ\';         end
if nargin < 5 || isempty(debug);     debug     = false;                 end
if nargin < 6 || isempty(poolcores); poolcores = 4;                     end

fhandles = {@maxtradepsec
    @sampleSpy
    @betacomponents
    @selrulecounts
    @rv
    @skewcomponents
    @countnullrets};

[hasFunc, pos] = ismember(fun, cellfun(@func2str,fhandles,'un',0));
if ~hasFunc
    error('Unrecognized function "%s".', fun)
end
fun             = fhandles{pos};
projectpath     = fileparts(mfilename('fullpath'));
[res, filename] = blockprocess(fun ,projectpath, varnames, cached,path2data,debug,poolcores, varargin{:});
end
%% Subfunctions

% Maximum trades per second
function res = maxtradepsec(s,cached)
% Count per id, date and second
Id            = cumsum([ones(1,'uint32'); diff(int8(s.data.Time(:,1))) < 0]);
Dates         = RunLength(s.mst.Date, double(s.mst.To-s.mst.From+1));
[un,~,subs]   = unique([Id, Dates, s.data.Time],'rows');
counts        = [un accumarray(subs,1)];
% Pick max per date
[date,~,subs] = unique(counts(:,2));
res           = table(date, accumarray(subs,counts(:,end),[],@max),'VariableNames',{'Date','Maxpsec'});
end

function res = betacomponents(s,cached, grid)
% DO NOT RELY on local id!

spdays = cached{2};
spyret = cached{1};
useon  = numel(cached) == 4;
if useon
    onret = cached{3};
end

% Dates and returns
dates = s.data.Datetime;
ret   = [NaN; diff(log(s.data.Price))];
idx   = [false; diff(rem(dates,1)) >= 0];
if useon
    [~,pos]   = ismembIdDate(s.mst.Permno, s.mst.Date, onret.Permno, onret.Date);
    ret(~idx) = onret.RetCO(pos);
else
    % Keep all except overnight
    ret = ret(idx,:);
end

if ~isempty(grid)
    % Mkt
    dt         = spyret{1}(:,1);
    [~,~,subs] = histcounts(mod(dt,1),grid);
    for ii = 1:numel(spyret)
        r          = spyret{ii}(:,2);
        spyret{ii} = accumarray(subs, nan2zero(r)+1,[],@(x) prod(x)-1);
    end

    % All others
    if ~issorted(s.mst.From)
        error('Unsorted mst!');
    end
    nmst       = size(s.mst,1);
    [subs,id]  = ndgrid(subs,1:nmst);
    [~,~,subs] = unique([id(:),subs(:)],'rows');
    ret        = accumarray(subs, nan2zero(ret)+1,[],@(x) prod(x)-1);
end

% Use a NaN when we don't have SPY returns
ngrid   = size(spyret{1},1);
spyret  = [NaN(ngrid,1); spyret];
days    = yyyymmdd2serial(double(s.mst.Date));
[~,pos] = ismember(days, spdays);
pos     = pos + 1;

% Map SP500 rets to stock rets
spret   = cat(1,spyret{pos});
prodret = spret.*ret;
subsID  = reshape(repmat(1:nmst,ngrid,1),[],1);
ikeep   = ~isnan(prodret);

% Store results
res     = s.mst(:,{'Permno','Date'});
res.Num = accumarray(subsID(ikeep), prodret(ikeep),[],[],NaN);
res.Den = accumarray(subsID(ikeep), spret(ikeep).^2,[],[],NaN);
end

function res = skewcomponents(s, cached)
res = s.mst(:,{'Date','Permno'});

% Add permno
nobs          = double(s.mst.To - s.mst.From + 1);
nmst          = size(s.mst,1);
subs          = RunLength((1:nmst)', nobs);
s.data.Permno = s.mst.Permno(subs);

% Dates and returns
s.data.Ret        = [NaN; diff(log(s.data.Price))];
ion               = [true; diff(fix(s.data.Datetime)) ~= 0] |...
                    [true; diff(s.data.Permno) ~= 0];
s.data.Ret(ion,1) = NaN;

% Use overnight
useon = numel(cached) == 2;
if useon
    s.data = addOvernightRet(s.data, cached{1});
end

% Filter out NaNs
ikeep = ~isnan(s.data.Ret);
subs  = subs(ikeep);
ret   = s.data.Ret(ikeep);

% Daily realized skewness: sqrt(nobs)*sum(r^3)/sum(r^2)^(3/2)
res.N   = accumarray(subs,1);
res.Sx3 = accumarray(subs, ret.^3);
res.Rv  = accumarray(subs, ret.^2);
end

function res = rv(s,cached)

% Dates and returns
dates = s.data.Datetime;
ret   = [NaN; diff(log(s.data.Price))];
idx   = [false; diff(rem(dates,1)) >= 0];

% Use overnight
useon = numel(cached) == 2;
if useon
    onret     = cached{1};
    [~,pos]   = ismembIdDate(s.mst.Permno, s.mst.Date, onret.Permno, onret.Date);
    ret(~idx) = onret.RetCO(pos);
else
    % Keep all except overnight
    ret = ret(idx,:);
end

% Number of observations per day
nobs = double(s.mst.To - s.mst.From + 1);
nmst = size(s.mst,1);
subs = RunLength((1:nmst)', nobs);

% Filter out NaNs
ikeep = ~isnan(ret);
subs  = subs(ikeep);
ret   = ret(ikeep);

% RV, sum and count
res    = s.mst(:,{'Permno','Date'});
res.RV = accumarray(subs, ret.^2,[],[],NaN);
res.Sx = accumarray(subs, ret   ,[],[],NaN);
res.N  = uint8(accumarray(subs,      1,[],[],NaN));
end

function res = countnullrets(s,cached)
res          = s.mst(:,{'Id','Permno','Date'});
idx          = mcolon(s.mst.From,1,s.mst.To);
subs         = RunLength(1:size(s.mst,1),s.mst.To-s.mst.From+1);
res.Nullrets = accumarray(subs(:), s.data.Price(idx),[],@(x) nnz(x(2:end)./x(1:end-1)==1));
end

%% OLD stuff
function res = selrulecounts(s,cached)
nfile  = uint16(cached{end});
vnames = {'Date','Val','Count'};

% Sort mst
if ~issorted(s.mst.From)
    s.mst = sortrows(s.mst,'From');
end
dates = RunLength(s.mst.Date, double(s.mst.To-s.mst.From+1));

% G127
[g127,~,subs] = unique([dates,s.data.G127_Correction(:,1)],'rows');
g127          = table(g127(:,1),g127(:,2),accumarray(subs,1), 'VariableNames', vnames);

% Correction
[correction,~,subs] = unique([dates, s.data.G127_Correction(:,2)],'rows');
correction          = table(correction(:,1),correction(:,2),accumarray(subs,1), 'VariableNames', vnames);

% Condition
[condition,~,subs] = unique(table(dates, s.data.Condition));
condition.Count    = accumarray(subs,1);
condition          = setVariableNames(condition,vnames);

% Null price
bins               = double(s.data.Price ~= 0);
[nullprice,~,subs] = unique([dates, bins],'rcachedmst.MedPrice== 0ows');
nullprice          = table(nullprice(:,1), nullprice(:,2), accumarray(subs,1),'VariableNames', vnames);

% Null size
bins              = double(s.data.Volume ~= 0);
[nullsize,~,subs] = unique([dates, bins],'rows');
nullsize          = table(nullsize(:,1), nullsize(:,2), accumarray(subs,1),'VariableNames', vnames);

res = {g127, correction, condition, nullprice, nullsize, nfile};
end