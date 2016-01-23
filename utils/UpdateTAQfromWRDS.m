addpath C:\Users\ok1011\Documents\github\wrds
wrds_install
w = wrds('olegkoma');
list = w.getDatasetNames('TAQ');

% Keep transactions only
idx = ~cellfun('isempty',regexp(list,'CT_|MAST_'));
list = list(idx);

% Download
topath = 'E:\TAQ\HFbetas\data\TAQ\raw';
for ii = 1:numel(list)
    filename = fullfile(topath, [list{ii}, '.zip']);
    libdtname = sprintf('TAQ.%s', list{ii});
    w.getDataset(libdtname, filename);
end