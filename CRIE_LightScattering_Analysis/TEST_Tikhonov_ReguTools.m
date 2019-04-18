clear all; close all

%Noisy data
[A,b_bar,x] = shaw(32);
randn('seed',41997);
e = 1e-3*randn(size(b_bar)); 
b = b_bar + e;

%CSVD decomposition
[U,s,V] = csvd(A);

%L-curve (Tikhonov regularization)
lambda_l = l_curve(U,s,b);   
axis([1e-3,1,1,1e3])
k_l = l_curve(U,s,b,'tsvd'); 
axis([1e-3,1,1,1e3])
