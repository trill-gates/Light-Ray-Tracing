function C = vdot(A,B),

% VDOT return dot product of rows of two [nx3] matrices. 
% i.e., an [nx3] matrix C with row n given by, dot(A(n,:), B(n,:))
%

C = A(:,1).*B(:,1) + A(:,2).*B(:,2) + A(:,3).*B(:,3);

return
