function [se,pval,sepw,pvalpw]=sharpeHAC(ret,type)
%This function performs HAC inference on the difference between 2 sharpe
%ratios. 
%Inputs:
    %ret= T*2 matrix of returns (type double).
    %type= (optional) specifies the kernel to be used to calculate Psi hat
    %2 options are available 'G' for the Parzen-Gallant kernel (default) and 'QS'
    %for the Quadratic Spectral kernel(type string).
%Outputs:
    %se= HAC standard error
    %pval= HAC p-value
    %sepw= HAC standard error pre-whitened
    %pvalpw= HAC p-value pre-whitened
    if not(ismember('type',who)), type='G'; end;
    srdiff = sharpeRatioDiff(ret);
    se = computeSE(ret,type);
    if nargout >=3
        sepw = computeSEpw(ret,type);
    end
    %calculating normal cdf recursively
    fun= @(x) (1/sqrt(pi*2))*exp(-0.5*x.^2);
    pval = 2 *integral(fun,-1000,-abs(srdiff)/se);
    if nargout == 4
        pvalpw = 2 * integral(fun,-1000,-abs(srdiff)/sepw);   
    end
end
