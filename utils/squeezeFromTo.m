function [stb, idx] = squeezeFromTo(tb)
% SQUEEZEFROMTO Consolidates contiguous From/To blocks
%
%   SQUEEZEFROMTO(TB) TB should be a table (or dataset) with 
%                     the fields 'From' and 'To'.

%                     If TB contains the field 'File', the 
%                     from/to indices are consolidated within 
%                     each file (grouping index).


if any(tb.From > tb.To) 
    error('squeezeFromTo:fromBiggerThanTo','''From'' cannot be bigger than ''To''.')
end
nrecs = size(tb,1);
if issorted(tb.From)
    if nargout == 2
        idx = 1:nrecs;
    end
else
    [tb, idx] = sortrows(tb,'From');
end

from   = [0; find(tb.To(1:end-1)+1 ~= tb.From(2:end))] + 1;
to     = [from(2:end) - 1; nrecs];

stb    = tb(from,:);
stb.To = tb.To(to);

end