function [Q,R] = givensqr(A)
% GIVENSQR - Given's Rotation factors m x n matrix A = QR where Q is m x m
%   orthogonal matrix (Q^{-1} = Q') and R is m x n upper triangular matrix.  
%
% Input: 
%   A = m x n matrix 
% 
% Output:
%   Q = m x m orthogonal matrix
%   R = m x n upper triangular matrix
% 
% See also houseqr, mgs, qr

% -------------------------------------------------------------------------
% Quinlan, J.
% 2022-12-31
% 
% -------------------------------------------------------------------------

    [m, n] = size(A);
    Q = eye(m);
    R = A;

    for j=1:n
        for i = m:-1:j+1
            G = eye(m);
            [c,s] = givens(R(i-1,j),R(i,j));
            G(i-1:i,i-1:i) = [c -s;s c];
            R = G'*R;
            Q = Q*G;
        end
    end
end

% Local Function
function [c,s] = givens(a,b)
    if abs(b) >= abs(a) 
        t = a/b;
        s = sign(b) / sqrt(1+t^2); 
        c = s*t;
    else
        t = b/a;
        c = sign(a) / sqrt(1 + t^2);
        s = c*t;
    end

end
