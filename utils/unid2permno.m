function permnos = unid2permno(unids)
% UNID2PERMNO Returns the permno corresponding to the UnID
%
%   UNID2PERMNO(UNIDS) Array of UnIDs
%
%
%   PERMNOS = UNID2PERMNO(...) Array of PERMNOS in preserved order.
%
% Note: only one permno per UnID (excluding NaN permnos). Also, some
%       UnID might be 65335, i.e. intmax('uint16').

% Preallocate
permnos      = NaN(size(unids));
% Load mapping
taq2crsp     = loadresults('taq2crsp');
taq2crsp     = unique(taq2crsp(:,{'permno','ID'}));
% Map
[idx,pos]    = ismember(unids, taq2crsp.ID);
permnos(idx) = taq2crsp.permno(pos(idx));
end