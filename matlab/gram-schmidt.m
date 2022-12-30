% Test Matrices 
% -------------------------------------------------------------------------
% A = [1 -1 4;1 4 -2;1 4 2;1 -1 0];
% A = [1 -2 -1;2 0 1;2 -4 2;4 0 0];
% A = [2    -3     1;  1    -2     1; 1    -3     2]
% A = [12 -51 4;6 167 -68;-4 24 -41];
A = [
  0.75126707 0.75126707 0.75126706 0.75126706;
  0.25509512 0.25509512 0.25509512 0.25509512;
  0.50595705 0.50595706 0.50595705 0.50595706;
  0.69907672 0.69907673 0.69907673 0.69907673;
  0.89090326 0.89090326 0.89090326 0.89090326
];


% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%% Gram-Schmitt Orthogonalization Process
%
% Given matrix A with n independent columns (e.g. rank(A) = n), return
% matrices Q and R with Q orthogonal and R upper triangular.  R contains
% inner products, i.e., R(i,j) = Q(:,i)'*A(:,j).
%
% Quinlan, J. 04/23/2022
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

clc

[m,n] = size(A);
Q = zeros(m,n);
R = zeros(n);
% P = zeros(m,n-1);
a = A(:,1);
 
R(1,1) = norm(a);
Q(:,1) = a ./ R(1,1);

for j = 2:n
    a = A(:,j);   
    % Compute projection
    for i = 1:j-1
        R(i,j) = Q(:,i)'*A(:,j);
        a = a - R(i,j)*Q(:,i);
    end
    R(j,j) = norm(a);
    Q(:,j) = a ./R(j,j);
end

% Print results
disp('A = ');disp(A);
disp('Q = ');disp(Q);
isothogonal = isortho(Q);
fprintf('Orthogonal = %d\n\n',isothogonal);
disp('R = ');disp(R);

QR = Q*R;
disp('QR = ');disp(QR);
disp('A = QR?');disp(abs(A - QR) < 2*eps);
% -------------------------------------------------------------------------





% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%% Modified Gram-Schmitt Orthogonalization Algorithm
%
% The above QR method does not generally produce accurate results in finite
% precision arithmetic (due to round-off when computing Q(:,i)).  
%
% A modified version gives better numerical accuracy by orthogonalizing
% A(:,2:n) to Q(:,1), thus modifiying the original set.  Then it computes 
% Q(:,2) and orthogonalizes A*(:,3:n) to Q(:,2) where A* is the modified
% vectors {a1^(1), a2^(1), ..., an^(1)}.
% 
% Quinlan, J. 04/23/2022
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


clc

[m,n] = size(A);
Q = zeros(m,n);
R = zeros(n);
a = A(:,1);

R(1,1) = norm(a);
Q(:,1) = a ./ R(1,1);

for j = 2:n
    a = A(:,j);   
    % Compute projection
    for i = 1:j-1
        R(i,j) = Q(:,i)'*a;
        a = a - R(i,j)*Q(:,i);
    end
    R(j,j) = norm(a);
    Q(:,j) = a ./R(j,j);
end

% Print results
disp('A = ');disp(A);
disp('Q = ');disp(Q);
isothogonal = isortho(Q);
fprintf('Orthogonal = %d\n\n',isothogonal);
disp('R = ');disp(R);

QR = Q*R;
disp('QR = ');disp(QR);

disp('A = QR?');disp(abs(A - QR) < eps);


%{
NOTES:
% -------------------------------------------------------------------------
Defintion: ||v|| = sqrt(<v,v>)  
     
Defintion: u and v are orthogonal is <u,v> = 0.  

Definition: An orthonormal set of vectors is an orthogonal set of unit
            vectors
     
Definition: Q is orthogonal if the columns form an orthonormal set
     
Theorem: An n x n matrix Q is orthongoal if and only if Q'Q=I


REFERENCES
---------------------------------------------------------------------------
  Golub, G. H., & Van Loan, C. F. (2012). Matrix computations (Vol. 3). JHU Press.
  Leon, S. (2011). Linear Algebra with Applications (8th).  Pearson.
%}
