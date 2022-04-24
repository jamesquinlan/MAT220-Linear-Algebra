%% Gram-Schmidt Orthogonalization Process
%
% Given matrix A with n independent columns (e.g. rank(A) = n), return
% matrices Q and R with Q orthogonal and R upper triangular.  R contains
% inner products, i.e., R(i,j) = Q(:,i)'*A(:,j).
%
% Quinlan, J. 04/23/2022

clc

%A = [1 -1 4;1 4 -2;1 4 2;1 -1 0];
A = [1 -2 -1;2 0 1;2 -4 2;4 0 0];
disp(A);

[m,n] = size(A);
Q = zeros(m,n);
R = zeros(n);
P = zeros(m,n-1);

% 
R(1,1) = norm(A(:,1));
Q(:,1) = A(:,1) ./ R(1,1);

for j = 2:n
       
    % Compute projection
    for i = 1:j-1
        R(i,j) = Q(:,i)'*A(:,j);
        P(:,j-1) = P(:,j-1) + R(i,j)*Q(:,i);
    end
    R(j,j) = norm(A(:,j) - P(:,j-1));
    Q(:,j) = (A(:,j) - P(:,j-1))./R(j,j);
end

disp(Q);
isortho(Q)
disp(R);




%% Modified Gram-Schmidt Orthogonalization Algorithm
%
% The above QR method does not generally produce accurate results in finite
% precision arithmetic (due to round-off when computing Q(:,i)).  
%
% A modified version gives better numerical accuracy by orthogonalizing
% A(:,2:n) to Q(:,1). 
% 
% Quinlan, J. 04/23/2022

clc
A = [1 -1 4;1 4 -2;1 4 2;1 -1 0];
% A = [1 -2 -1;2 0 1;2 -4 2;4 0 0];
[m,n] = size(A);
Q = zeros(m,n);
R = zeros(n);
AA = A;  % NOTE: A changes within loop so AA is a copy

for i = 1:n
    R(i,i) = norm(A(:,i));
    Q(:,i) = A(:,i) / R(i,i);
    for j = i+1:n
       R(i,j) = Q(:,i)'*A(:,j);
       A(:,j) = A(:,j) - R(i,j)*Q(:,i);
    end
end

% Print results
disp('A = ');disp(A);
disp('Q = ');disp(Q);
isothogonal = isortho(Q);
fprintf('Orthogonal = %d\n\n',isothogonal);
disp('R = ');disp(R);

QR = Q*R;
disp('A = QR?');disp(abs(AA - QR) < eps);


%{
NOTES:
-------------------------------
Defintion: ||v|| = sqrt(<v,v>)  
     
Defintion: u and v are orthogonal is <u,v> = 0.  

Definition: An orthonormal set of vectors is an orthogonal set of unit
            vectors
     
Definition: Q is orthogonal if the columns form an orthonormal set
     
Theorem: An n x n matrix Q is orthongoal if and only if Q'Q=I


REFERENCES:
-------------------------------
  Golub, G. H., & Van Loan, C. F. (2012). Matrix computations (Vol. 3). JHU Press.
  Leon, S. (2011). Linear Algebra with Applications (8th).  Pearson.
%}

