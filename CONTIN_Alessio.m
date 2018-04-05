% CREATED: Alessio Caciagli, University of Cambridge, March 2018
clear variables;close all
%% Input

Ng = 161; %Number of grid points
if bitget(Ng,1) == 0 %(even)
    Ng = Ng + 1;
end

%Fictitious data
t = logspace(-7,2,100);
y = exp(-t/1e-3);

%Grid points
s_min = 1/t(end);
s_max = 1/t(1);

% Frequency spectrum: f_min = 1/t_max and viceversa
% Samplings
%if (isempty(varargin{1}) || strcmp(varargin{1},'logarithmic'))
    s = linspace(log10(s_min),log10(s_max),Ng);
%elseif strcmp(varargin{1},'linear')
%    s = linspace(s_min,s_max,Ng);
%else
%    error('Unknown sampling.');
%end
c = [4*ones(1,(Ng-3)/2);2*ones(1,(Ng-3)/2)]; %Simpson rule for coeff
c = c(:);
h = (s(end) - s(1))/(Ng-1);
c = (h/3)*[1;c;4;2];
%% Input post-process

t = t(:); % must be column
y = y(:); % must be column
s = s(:); % must be column
c = c(:); % must be column
%g0 = g0(:); % must be column

[sM,~] = meshgrid(s,t);
[cM,tM] = meshgrid(c,t);
A = cM.*(log(10)*10.^(sM)).*exp(-tM.*10.^(sM));


c = [ones(1,(Ng-3)/2);4*ones(1,(Ng-3)/2)]; %Simpson rule for coeff


