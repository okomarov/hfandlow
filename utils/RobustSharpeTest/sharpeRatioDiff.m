function srdiff = sharpeRatioDiff(ret)
% Computes the difference betweeen two Sharpe ratios
% Inputs:
    % ret = T*2 matrix of returns (type double)
% Outputs:
    % diff = difference of the two Sharpe ratios
% Note:
    % returns are assumed to be in excess of the risk-free rate already
    
    SR = mean(ret,1)./std(ret,[],1);
    srdiff = SR(1) - SR(2);
end
