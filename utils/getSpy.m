function spy = getSpy(freq, from, to)
% spy = getSpy(freq)
%
% Note: SPY has permno 84398

if nargin < 2 || isempty(from), from = 0;   end
if nargin < 3 || isempty(to),   to   = inf; end

% Sampling params
if nargin < 1 || isempty(freq), freq = 5; end
name = sprintf('spy%dm.mat',freq);

% Unsampled
if isinf(freq)
    spy = getTaqData('symbol','SPY',from,to);

% Sampled
else
    try
        spy = loadresults(name);
    catch
        % Sample
        fprintf('%s: sampling SPY at %d min.\n', mfilename, freq)
        opt.grid = (9.5/24:freq/(60*24):16/24)';
        writeto  = '.\results';
        
        % Sample at x min
        [spy, filename] = Analyze('sampleSpy',[],[],[],[],opt);
        
        % Rename to append the sampling frequency
        newfullname = fullfile(writeto, name);
        movefile(fullfile(writeto,filename), newfullname);
    end
    % Filter dates
    idx = in(serial2yyyymmdd(spy.Datetime),[from, to]);
    spy = spy(idx,:);
end
end