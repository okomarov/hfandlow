% Get sp100 (Gvkeyx == 664) constitutents

% Get LC, LU, LN linktype
ccm           = loadresults('ccm');
ccm           = ccm(ismember(ccm.Linktype,[1,2,6]),:);
ccm           = ccm(ccm.Linkprim ~= 'J',:);
ccm           = sortrows(ccm,{'Lpermno','Linkdt'});
idx           = find(isfeatchange(ccm(:,[5,1,7])));
enddate       = ccm.Linkenddt([idx(2:end)-1; min(idx(end), size(ccm,1))]);
ccm           = ccm(idx,:);
ccm.Linkenddt = enddate;

const = loadresults('constituents');
const = const.constit(:,1:4);
const = const(const.Gvkeyx == 664,:);
const = unique(const);

% First XS intersection
idx = ismember(ccm.Gvkey,const.Gvkey);
ccm = sortrows(ccm(idx,:),'Gvkey');

% TS intersection
from  = 19930101;
to    = serial2yyyymmdd(now);
const = const(const.Thru >= from,:);
ccm   = ccm(ccm.Linkenddt >= from,:);
const = const(const.From <= to,:);
ccm   = ccm(ccm.Linkdt <= to,:);

const.From(const.From < from)     = from;
ccm.Linkdt(ccm.Linkdt < from)     = from;
const.Thru(const.Thru > to)       = to;
ccm.Linkenddt(ccm.Linkenddt > to) = to;

% Unstack
const = pivotFromTo(const(:,[1,3,4,2]));
ccm   = pivotFromTo(ccm(:,[1,7,8,5]));

% Ensure same XS
[idx,pos]   = ismember(const.Id, ccm.Id);
const.Panel = const.Panel(:, idx);
ccm.Panel   = ccm.Panel(:, pos(idx));
ccm.Id      = ccm.Id(pos(idx));

permnos                            = ccm.Panel{:,:};
permnos(logical(const.Panel{:,:})) = 0;

% Stack back
[idate,~]  = find(permnos);
sp100const = table(ccm.Date(idate), permnos(permnos~=0),...
                   'VariableNames',{'Date','Permno'});
% Save
filename = sprintf('%s_sp100const.mat',datestr(now,'yyyymmdd_HHMM'));
save(fullfile('results', filename), 'sp100const')
     