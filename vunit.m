function n = vunit(n)

% VUNIT Return a matrix so that the rows are normalized to 
% unit length i.e. each row has norm = 1.
%

% normalization factor
n0 = n.^2;
n0 = n0';
n0 = sum(n0);
n0 = sqrt(n0);
n0 = n0';

% divide each column
n = n./ n0(:,ones([1 size(n,2) ])); 
