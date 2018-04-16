function Q = CONTIN_Alessio_Q(A)
%CONTIN_Alessio_Q Computes the Q matrix of the factorization QR of matrix A
%
% CONTIN_Alessio_Q(A)
%
% Given a matrix A (m x n, with m > n), CONTIN_Alessio_qr computes the Q
% matrix resulting from the QR factorization of A.

[m,n] = size(A);
if m < n
    error('Error in QR factorization: the matrix must have more rows than columns!')
end
U = zeros(m,n);
Q = eye(m);
for i = 1:n
    [~,beta,u] = gen_hh(A(i:end,i));
    A = app_hh(A,beta,[zeros(i-1,1);u]);
    U(i:m,i)=sqrt(beta)*u;
end
for i = n:-1:1
    Q = app_hh(Q,1,U(:,i));
end

