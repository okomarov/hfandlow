function [betas, permnos] = getBetas(freq, useon, useproxy, issp, iscs, keeplong)
% GETBETAS Loads betas, applies sp500 and/or common shares filers, unstacks
%
%   getBetas(lookback, freq, useon, useproxy, issp, iscs, keeplong)
%
%   NOTE: Does NOT fill in between NaNs

if nargin < 1 || isempty(freq),     freq     = 5;       end
if nargin < 2 || isempty(useon),    useon    = true;    end
if nargin < 3 || isempty(useproxy), useproxy = false;   end
if nargin < 4 || isempty(issp),     issp     = true;    end
if nargin < 5 || isempty(iscs),     iscs     = true;    end
if nargin < 6 || isempty(keeplong), keeplong = false;   end

% Load all betas
betas = estimateBetaComponents(freq, useon, useproxy);

% Filter for sp500 members
if issp
    fprintf('%s: filtering for sp500 members.\n', mfilename)
    idx   = issp500member(betas(:,{'Permno','Date'}));
    betas = betas(idx,:);
end

% Filter for common share (share type code 10 and 11)
if iscs
    fprintf('%s: filtering for common shares only.\n', mfilename)
    idx   = iscommonshare(betas(:,{'Permno','Date'}));
    betas = betas(idx,:);
end

if nargout == 2
    permnos = unique(betas.Permno);
end

% Unstack betas
if ~keeplong
    betas = unstack(betas(:,{'Date','Permno','Beta'}), 'Beta','Permno');
    betas = sortrows(betas,'Date');
end

% Convert to double
betas = tbextend.varfun(@double, betas, 'RenameVariables', false);
end