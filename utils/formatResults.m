function out = formatResults(coeff, se, pval, colheaders, rowheaders)
% formatResults(coeff, se, pval, colheaders, rowheaders)

% Optional arguments
coeffFmt   = '%.2f';
seFmt      = '{(%.3f)}';
seFontSize = '\footnotesize';
txtBefore  = ['\begin{tabu}{l *{#nchead#}{S[table-format   = -2.2,', 10,...
              '                      table-space-text-post = $^{****}$,',10,...
              '                      table-align-text-post = false,',10,...
              '                      table-text-alignment  = center]',10,...
              '                      @{}',10,...
              '                      }',10,...
              ' }',10,...
              '\toprule'];
txtAfter   = ['\bottomrule',10,...
              '\end{tabu}'];

          
% Pre-allocate body
sz   = size(coeff);
body = cell(sz(1)*2, sz(2));

% Format coefficients and se
cellcoeff       = arrayfun(@(x)sprintf(coeffFmt,x), coeff, 'un', 0);
body(2:2:end,:) = arrayfun(@(x)sprintf(seFmt,x), se, 'un', 0);

% Add *, **, *** for significance at 10, 5 and 1%
body(1:2:end,:) = addStars(cellcoeff, pval);

% Clean out NaNs
inan       = logical(kron(isnan(coeff),ones(2,1)));
body(inan) = {[]};

% improve readability
body = fixedWidthColSpacing(body);

% Column headers and midrule
nchead          = numel(colheaders);
[chead,midrule] = deal(cell(1,sz(2)));
colperhead      = sz(2)/nchead;
for ii = 1:nchead
    chead{ii}   = sprintf('\\multicolumn{%d}{c}{%s}',colperhead, colheaders{ii});
    midrule{ii} = sprintf('\\cmidrule(lr){%d-%d}',(ii-1)*colperhead+2, ii*colperhead+1);
end
chead{nchead} = strrep(chead{nchead},'&','\\');

% Add row headers
rhead            = [rowheaders; repmat({sprintf('\\rowfont{%s}',seFontSize)},1,sz(1))];
rhead            = [cell(2,1); rhead(:)];
% Improve readability
rhead([1,3:end]) = fixedWidthColSpacing(rhead([1,3:end]));

% Add endlines
body([1:2:end,end],end) = padEntry(body([1:2:end,end],end), ' \\');
body(2:2:end-1,end)     = padEntry(body(2:2:end-1,end), ' \\[3pt]');
chead(end)              = padEntry(chead(end), ' \\');

% Combine all together
out = [rhead [chead; midrule; body]];

% Add column separators
out(1,1:nchead)    = padEntry(out(1,1:nchead), ' &');
out(3:end,1:end-1) = padEntry(out(3:end,1:end-1), ' &' );

% Concatenate into one string
out(:,end) = padEntry(out(:,end), char(10));
out        = num2cell(out,1);
out        = strcat(out{:});
out        = [strrep(txtBefore,'#nchead#',sprintf('%d',sz(2))),...
              [out{:}],...
              txtAfter];
end

function c = addStars(c, pval)
a      = [-inf, 0.01, 0.05, 0.1];
starts = {'$^{***}$','$^{**}$','$^*$'};
for ii = 1:3
    idx    = a(ii) < pval & pval <= a(ii+1);
    c(idx) = strcat(c(idx), starts(ii));
end
end

function c = padEntry(c, padval)
for ii = 1:numel(c)
    c{ii} = [c{ii}, padval];
end
end

function c = fixedWidthColSpacing(c)
% Add fixed width separation
len    = cellfun(@length, c);
maxlen = max(len(:));
for ii = 1:numel(c)
    c{ii} = sprintf('%s%s',repmat(' ',1, maxlen-len(ii)), c{ii});
end
end