function [R,c] = CONTIN_Alessio_qr(A,b)
%CONTIN_Alessio_qr Computes the qr factorization of the problem Ax=b (x
%unknown)
%
% CONTIN_Alessio_qr(A,b)
%
% Given a matrix A (m x n, with m > n) and a vector b (m x 1), 
% CONTIN_Alessio_qr computes the QR factorization of the problem Ax = b (x 
% unknown) by generating the matrix R and vector c so that they fulfill
% Rx=c, with R (n x n) upper triangular matrix.

[m,n] = size(A);
if m < n
    error('Error in QR factorization: the matrix must have more rows than columns!')
end

for i = 1:n
    [~,beta,u] = gen_hh(A(i:end,i));
    A = app_hh(A,beta,[zeros(i-1,1);u]);
    b = app_hh(b,beta,[zeros(i-1,1);u]);
end
R = A(1:n,1:n);
c = b(1:n);