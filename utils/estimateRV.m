function rv = estimateRV(lookback, freq, useon)
if nargin < 1 || isempty(lookback), lookback = 1;     end
if nargin < 2 || isempty(freq),     freq     = 5;     end
if nargin < 3 || isempty(useon),    useon    = true; end

writeto = '.\results\';

try
    name = matname('rv',freq, useon);
    rv   = loadresults(name);
catch
    
    % Sample if data doesn't exist
    path2data = sprintf('.\\data\\TAQ\\sampled\\%dmin', freq);

    % Add overnight return or discard
    if useon
        load(fullfile(path2data,'master'),'-mat','mst');
        reton                 = loadresults('return_intraday_overnight');
        [idx,pos]             = ismembIdDate(reton.Permno,reton.Date, mst.Permno, mst.Date);
        mst.RetCO(pos(idx),1) = reton.RetCO(idx);
        fun                   = @(x) {mst(x,{'Permno','Date','RetCO'})};
        mst                   = accumarray(mst.File,(1:size(mst))',[], fun);
    end
       
    % Calculate beta components: sum(r*benchr) and sum(benchr^2)
    fprintf('%s: creating RV components at %d min.\n', mfilename, freq)
    [rv,filename] = Analyze('rv', [], mst, path2data);
    
    % Rename to append the sampling frequency
    name        = regexp(filename,'\w+?(?=\.mat)','match','once');
    name        = [matname(name,freq, useon),'.mat'];
    newfullname = fullfile(writeto, name);
    movefile(fullfile(writeto,filename), newfullname);
end
fprintf('%s: calculating RV with %d day lookback.\n', mfilename, lookback)

% Sort (composite key for speed) and create subs
key        = uint64(rv.Permno) * 1e8 + uint64(rv.Date);
[~,isrt]   = sort(key);
rv         = rv(isrt,:);
[~,~,subs] = unique(rv.Permno);

% RV = ?x^2; 
% RVscaled = RV/n = E[x^2]; 
% Var = E[x^2] - E[x]^2 = RVscaled - (?x/n)^2;
RVroll      = accumarray(subs,        rv.RV, [], @(x) runsum(lookback,x));
Sx          = accumarray(subs,        rv.Sx, [], @(x) runsum(lookback,x));
N           = accumarray(subs, double(rv.N),  [], @(x) runsum(lookback,x));
rv.RVscaled = cat(1,RVroll{:})./cat(1,N{:});
rv.Var      = rv.RVscaled - (cat(1,Sx{:})./cat(1,N{:})).^2;
end

function name = matname(name, freq, useon)
if useon, useon = 'on'; else useon = ''; end
name = sprintf('%s%dm%s%s', name, freq, useon);
end