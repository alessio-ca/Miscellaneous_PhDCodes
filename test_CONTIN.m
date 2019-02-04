% Target g(s) function:
g0 = [0 0 10.1625 25.1974 21.8711 1.6377 7.3895 8.736 1.4256 0 0]';
% s space:
s0  = logspace(-3,6,length(g0))';

semilogx(s0,g0); title('Target times distribution - Press any key...'); shg
pause

t = linspace(0.01,500,100)';
[sM,tM] = meshgrid(s0,t);
A = exp(-tM./sM);
% Data:
data = A*g0 + 0.07*rand(size(t));

alpha = 1e-2;

% Guess s space and g function
s = logspace(-3,6,20)';
g = ones(size(s));
[g,yfit,cfg] = rilt(t,data,s,g,alpha,'logarithmic',[],[],[],{'g>0'},[],[]);

subplot(1,2,1)
plot(t,data,'.',t,yfit,'ro'); title('data and fitting')
subplot(1,2,2)
semilogx(s0,g0/max(g0),s,g/max(g),'o-'); title('g-target and g');