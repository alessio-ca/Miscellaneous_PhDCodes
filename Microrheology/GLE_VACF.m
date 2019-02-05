function [vacf,k,err] = GLE_VACF(tau,omega,lambda,gamma_l,varargin)

% [vacf,k,err] = GLE_VACF(tau,omega,lambda,gamma_l) Evaluates the analytical 
% VACF of a particle in a viscoelastic fluid subjected to an harmonic potential, 
% according to the GLE equation. 
% The memory kernel of the viscoelastic medium follows a power law (see [1]). 
%
%   TAU is a vector of time values.
%   OMEGA is the pulsation of the harmonic potential
%   LAMBDA and GAMMA_L parametrize the memory kernel (See [1]).
%   K is the maximum term considered in the Mittag-Leffler expansion.
%   ERR is the absolute error on the result by truncation of the series.
%
%  [vacf,k,err] = GLE_VACF(tau,omega,lambda,gamma_l,'PropertyName',PropertyValue) permits
%   to set the value of PropertyName to PropertyValue.
%   Admissible Properties are:
%       k_max       -   Maximum term in the power expansion   

% REFERENCES
%
%   [1] Desposito and Vinales, Physical Review E - Statistical, Nonlinear,
%   and Soft Matter Physics, 2009, 80, 2, 1-7

% Maximum term in pow exp
k_max = +Inf;
for n = 1:2:length(varargin)
    if strcmpi(varargin{n},'k_max')
        k_max = varargin{n+1};
    end
end

%Single term in the VACF series
func = @(x,k) (-1).^k * ml_func([2-lambda,1+2*k,k+1,1,((omega*x).^(2*k))],-gamma_l*x^(2-lambda));

k = zeros(size(tau));
vacf = zeros(size(tau));
err = ones(size(tau));
for i=1:length(tau)
    %Truncation method for the series evaluation
    while (err (i)> eps && k(i) < k_max)
        temp = func(tau(i),k(i));
        vacf(i) = vacf(i) + temp;
        k(i) = k(i) + 1;
        err(i) = abs(temp);
    end
end