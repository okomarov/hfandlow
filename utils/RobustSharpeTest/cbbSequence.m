function [sequence] = cbbSequence(T,b)
% Computes a circular block bootstrap sequence applied to (1:T)
% Inputs: 
    % T
    % b = block size
% Outputs:
    % sequence = bootstrap sequence
% Note:
    % 
    
    l             = floor(T/b)+1;
    indexSequence = [1:T,1:b]';
    startPoints   = randi(T,1,l);
    pos           = bsxfun(@plus,startPoints,(0:b-1)');
    sequence      = indexSequence(pos(:));
    sequence      = sequence(1:T);
end
