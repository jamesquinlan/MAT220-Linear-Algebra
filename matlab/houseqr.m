function [Q,R] = houseqr(A)
% HOUSEQR - Householder Reflection factors m x n matrix A = QR where Q is m x m
%   orthogonal matrix (Q^{-1} = Q') and R is m x n upper triangular matrix.
%
% Input:
%   A = m x n matrix
%
% Output:
%   Q = m x m orthogonal matrix
%   R = m x n upper triangular matrix
%
% See also givensqr, mgs, qr

% -------------------------------------------------------------------------
% Quinlan, J.
% 2022-12-31
%
% -------------------------------------------------------------------------

[m,n]=size(A);
Q=eye(m);
R=A;
for j=1:n
    v = house(R(j:m,j));
    c = 2/(v'*v);
    R(j:m,j:n) = R(j:m,j:n)-c*v*(v'*R(j:m,j:n));
    Q(:,j:m) = Q(:,j:m)-c*(Q(:,j:m)*v)*v';
end

function v=house(x)
v=x;
v(1)=v(1)+sign(v(1))*norm(x);


%{
A = [
0.751267 0.751267 0.751267 0.751267 ;
0.255095 0.255095 0.255095 0.255095 ;
0.505957 0.505957 0.505957 0.505957 ;
0.699077 0.699077 0.699077 0.699077 ;
0.890903 0.890903 0.890903 0.890903]

%}
