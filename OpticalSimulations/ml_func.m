function f=ml_func(aa,z,varargin)
% f = ML_FUNC(aa,z) Evaluates the 5-parameter Mittag-Leffler function of Z
%   and its derivatives.
%   AA is a vector of coefficients (at least one must be specified).
%
% f = ML_FUNC(aa,z,'PropertyName',PropertyValue) permits to set the value 
%   of PropertyName to PropertyValue.
%   Admissible Properties are:
%       n       -   Order of the derivative  

% NOTE
% The 5-parameter Mittag-Leffler function E_{a,b,gamma,q,p} (my invention!) 
% is defined as:
% E = sum_{k=0}^{infty} p*Gamma(gamma+k*q)*z^k/Gamma(gamma)/k!/Gamma(alpha*k+beta)

%Parameter assignation
aa = [aa,1,1,1,1];
a=aa(1);
b=aa(2);
c=aa(3);
q=aa(4);
p=aa(5);

f=0;
k=0;
fa=1;
% Deriv order
n = 0;
for i = 1:2:length(varargin)
    if strcmpi(varargin{i},'n')
        n = varargin{i+1};
    end
end

%Precision target
eps0=eps;

if n==0
    %Truncation method for the series evaluation (ML series)
    while norm(fa,1) >= eps0
        fa = (p*gammaAcc(k*q+c)*z.^k)./(gammaAcc(a*k + b)*gammaAcc(c)*factorial(k));
        f=f+fa;
        k=k+1;
    end
    %If convergence is not reached, throw error
    if any(~isfinite(f))
        eps1=round(-log10(eps0));
        % For the 2-parameter ML, attempt the use of a slower (but
        % accurate)formula [1].
        if c*q==1
            f=mlf(a,b,z,eps1);
            f=reshape(f,size(z));
        else
            error('Error: truncation method failed')
        end
    end
else
    %Recursive formula for the derivative evaluation
    aa(2)=aa(2)+n*aa(1);
    aa(3)=aa(3)+aa(4)*n;
    f = p*(gammaAcc(q*n+c)/gammaAcc(c))*ml_func(aa,z,'n',0);
end
    

