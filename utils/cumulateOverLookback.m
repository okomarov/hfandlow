function betas = cumulateOverLookback(betas,lookback)
if ischar(lookback)

    switch lookback
        case {'month','m'}
            [~,pos,subs] = unique(uint64(betas.Permno)*1e6 + uint64(betas.Date/100));
            out          = betas(pos,{'Permno','Date'});
            tmp          = accumarray(subs, betas.Num,[],@nansum)./ accumarray(subs, betas.Den,[],@nansum);
            out.Beta     = tmp;

        case {'year','y'}
            myunstack = @(tb,vname) sortrows(unstack(betas(:,{'Date','Permno',vname}),vname,'Permno'),'Date');
            num       = myunstack(betas,'Num');
            den       = myunstack(betas,'Den');
            dates     = num.Date;
            permnos   = num.Properties.VariableNames(2:end);
            num       = num{:,2:end};
            den       = den{:,2:end};

            [~,pos] = unique(dates/100,'last');

            num = cumsum(nan2zero(num));
            num = num(pos,:);
            den = cumsum(nan2zero(den));
            den = den(pos,:);

            pad = NaN(11,numel(permnos));
            num = [pad; num(13:end,:) - num(1:end-12,:)];
            den = [pad; den(13:end,:) - den(1:end-12,:)];
            
            out = num./den;
    end
    
elseif lookback == 1
    betas.Beta = betas.Num./betas.Den;

else
    fprintf('%s: calculating betas with %d day lookback.\n', mfilename, lookback)

    % Sort (composite key for speed) and create subs
    key        = uint64(betas.Permno) * 1e8 + uint64(betas.Date);
    [~,isrt]   = sort(key);
    betas      = betas(isrt,:);
    [~,~,subs] = unique(betas.Permno);

    % Beta
    num        = accumarray(subs, betas.Num, [], @(x) runsum(lookback,x));
    den        = accumarray(subs, betas.Den, [], @(x) runsum(lookback,x));
    betas.Beta = cat(1,num{:})./cat(1,den{:});
end
end