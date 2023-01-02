function x = backsub(A,b)
% BACKSUB backward substitution algorithm. 
%
% INPUT:
%   A = upper triangular matrix
%   b = vector
% 
% SYNOPSIS (Given the system of equations)
% a_11 x_1 + a_12 x_2 + ... + a_{1,n-1} x_{n-1} + a_{1,n} x_n = b_1
%            a_22 x_2 + ... + a_{2,n-1} x_{n-1} + a_{2,n} x_n = b_2
%                   .                                         .
%                    .                                        .
%                     .                                       .
%                         a_{n-1,n-1} x_{n-1} + a_{n-1,n}x_n = b_{n-1}
%                                                 a_{n,n}x_n = b_{n}
%
% See also forsub

% Error checking
if ~all(triu(A')-diag(diag(A)) ==0,'all')
    error('Matrix must be triangular. First try Gaussian Elimination');
end
if any(diag(A)==0)
    error('Matrix is singular'); 
end

n = size(A,2);
x = zeros(n,1);


% Backsub process
for i = n:-1:1
    x(i) = (1/A(i,i))*(b(i) - sum(A(i,i+1:n)*x(i+1:n)));
end
    
end

%% Notes
%{
    If A(i,i) is small, then 1/A(i,i) is large and will magnify error.
    Can we precondition matrix so to eliminate this?  
    Consider the example below and notice the difference in vector values.

Example:
------------
A=rand(50);
A=triu(A);
b=rand(50,1);
x1=A\b;
x2=backsub(A,b);
abs(x1(1) - x2(1)); % or 
abs(x1 - x2);  % how to justify this discrepency
------------

%}