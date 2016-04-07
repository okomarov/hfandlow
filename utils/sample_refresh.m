function [iA, iS] = sample_refresh(A,S)
% date = 20050722;
% A = getTaqData('symbol','AAPL',date,date);
% S = getSpy(inf,date,date);
% A = unique(A.Datetime);
% S = unique(S.Datetime);

nA = numel(A);
nS = numel(S);
iA = false(nA,1);
iS = false(nS,1);
cA = 1;
cS = 1;

while cA < nA && cS < nS
    dateA = A(cA);
    dateS = S(cS);

    if dateA > dateS
        iA(cA) = true;
        while S(cS) <= dateA
            cS = cS + 1;
        end
        iS(cS-1) = true;
        
        cA = cA+1;
        
    else
        iS(cS) = true;
        while A(cA) <= dateS
            cA = cA + 1;
        end
        iA(cA-1) = true;
        
        cS = cS+1;
    end
end
end
