function signal_filtered = filter_signal(mdl, signal, useMatlab)
if nargin < 3
    useMatlab = true;
end

[nrow,nser,nsig] = size(signal);
signal_filtered  = NaN(nrow,nser,nsig);
PQ               = int8([mdl.P, mdl.Q]);
for ii = 1:nsig
    ivalid = ~isnan(signal(:,:,ii));
    parfor c = 1:nser
        idx = ivalid(:,c);
        if nnz(idx) <= 10
            continue
        end

        out = NaN(nrow,1);
        y   = signal(idx,c,ii);

        if useMatlab
            fit    = mdl.estimate(y,'Display','off');
            fitted = y - fit.infer(y);
        else
            fitted = fitWithPython(y, PQ, mdl);
        end

        out(idx)                = fitted;
        signal_filtered(:,c,ii) = out;
    end
end
end

function yhat = fitWithPython(x, PQ, mdl)
x2 = py.numpy.array(x');
P  = double(PQ(1));
try
    fit  = py.statsmodels.tsa.arima_model.ARMA(x2,PQ).fit(pyargs('method','mle'));
    res  = py.statsmodels.tools.tools.maybe_unwrap_results(fit);
    yhat = double(py.array.array('d',res.fittedvalues))';
    return
catch
end

try
    xprev  = x(1:end-P);
    params = [x(1+P:end)\[ones(numel(xprev),1), xprev], 0];
    fit    = py.statsmodels.tsa.arima_model.ARMA(x2,PQ).fit(pyargs('method','mle','start_params',params));
    res    = py.statsmodels.tools.tools.maybe_unwrap_results(fit);
    yhat   = double(py.array.array('d',res.fittedvalues))';
catch
    fit  = mdl.estimate(x,'Display','off');
    yhat = x - fit.infer(x);
end
end
