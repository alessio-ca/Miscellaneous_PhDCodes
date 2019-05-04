function [reg_c,rho_c,eta_c] = CRIE_corner(rho,eta,reg_param,M,N)
% CRIE_Corner Locate the "corner" of the L-curve from CRIE analysis.
%
% [reg_c,rho_c,eta_c] =
%        CRIE_corner(rho,eta,reg_param)
%        CRIE_corner(rho,eta,reg_param,M)
%
% Locates the "corner" of the L-curve in log-log scale.
%
% It is assumed that corresponding values of || A x - b ||, || L x ||,
% and the regularization parameter are stored in the arrays rho, eta,
% and reg_param, respectively.
%
% An fourth argument M specifies an upper bound for eta, below which
% the corner should be found.

% Alessio Caciagli, University of Cambridge, April 2019.

% Ensure that rho and eta are column vectors.
rho = rho(:); eta = eta(:); reg_param = reg_param(:);

% Order rho in reverse order
if (rho(1) < rho(end))
    rho = flipud(rho);
    eta = flipud(eta);
    reg_param = flipud(reg_param);
end
    
    

% MODIFICATION: ensure constant sampling density (logspace for rho)
logSample = logspace(log10(rho(1)),log10(rho(end)),length(rho))';
rhoSample = zeros(size(rho));
etaSample = zeros(size(rho));

for i=1:length(logSample)
    [~,rpi] = min((rho-logSample(i)).^2);
    rhoSample(i) = rho(rpi);
    etaSample(i) = eta(rpi);
end
rho= unique(rhoSample);
eta = flipud(unique(etaSample));

% Set default reg_param vector method.
if (nargin==2), reg_param = (1:length(rho))'; end

% Set this logical variable to 1 (true) if the corner algorithm
% should always be used, even if the Spline Toolbox is available.
alwayscorner = 0;

% Set default parameters for treatment of discrete L-curve.
deg   = 2;  % Degree of local smooting polynomial.
q     = 2;  % Half-width of local smoothing interval.
order = 4;  % Order of fitting 2-D spline curve.

% Initialization.
if (length(rho) < order)
  error('Too few data points for L-curve analysis')
end

% Restrict the analysis of the L-curve according to M (if specified).
if (nargin>=4)
    index = find(rho < M);
    if ~isempty(index)
        rho = rho(index); eta = eta(index); reg_param = reg_param(index);
    end
end

% Restrict the analysis of the L-curve according to N (if specified).
if (nargin==5)
    index = find(eta < N);
    if ~isempty(index)
        rho = rho(index); eta = eta(index); reg_param = reg_param(index);
    end
end

% Use the adaptive pruning algorithm to find the corner, if the
% Spline Toolbox is not available.
if ~exist('splines','dir') || alwayscorner
    %error('The Spline Toolbox in not available so l_corner cannot be used')
    reg_c = corner(rho,eta);
    rho_c = rho(reg_c);
    eta_c = eta(reg_c);
    return
end

% Otherwise use local smoothing followed by fitting a 2-D spline curve
% to the smoothed discrete L-curve.

% Convert to logarithms.
lr = length(rho);
lrho = log(rho); leta = log(eta); slrho = lrho; sleta = leta;

% MODIFICATION: Scale data to have same boundaries
% scalerho = abs(lrho(end) - lrho(1));
% scaleeta = abs(leta(end) - leta(1));
% 
% scale = scalerho/scaleeta;
% if scale > 1
%     leta = leta*scale;
%     sleta = sleta*scale;
% else
%     lrho = lrho*(1/scale);
%     slrho = slrho*(1/scale);
% end



% For all interior points k = q+1:length(rho)-q-1 on the discrete
% L-curve, perform local smoothing with a polynomial of degree deg
% to the points k-q:k+q.
v = (-q:q)'; A = zeros(2*q+1,deg+1); A(:,1) = ones(length(v),1);
for j = 2:deg+1, A(:,j) = A(:,j-1).*v; end
for k = q+1:lr-q-1
    cr = A\lrho(k+v); slrho(k) = cr(1);
    ce = A\leta(k+v); sleta(k) = ce(1);
end

% Fit a 2-D spline curve to the smoothed discrete L-curve.
sp = spmak((1:lr+order),[slrho';sleta']);
pp = ppbrk(sp2pp(sp),[4,lr+1]);

% Extract abscissa and ordinate splines and differentiate them.
% Compute as many function values as default in spleval.
P     = spleval(pp);  dpp   = fnder(pp);
D     = spleval(dpp); ddpp  = fnder(pp,2);
DD    = spleval(ddpp);
ppx   = P(1,:);       ppy   = P(2,:);
dppx  = D(1,:);       dppy  = D(2,:);
ddppx = DD(1,:);      ddppy = DD(2,:);

% Compute the corner of the discretized .spline curve via max. curvature.
% No need to refine this corner, since the final regularization
% parameter is discrete anyway.
% Define curvature = 0 where both dppx and dppy are zero.
k1    = dppx.*ddppy - ddppx.*dppy;
k2    = (dppx.^2 + dppy.^2).^(1.5);
I_nz  = find(k2 ~= 0);
kappa = zeros(1,length(dppx));
kappa(I_nz) = -k1(I_nz)./k2(I_nz);
[kmax,ikmax] = max(kappa);
x_corner = ppx(ikmax); y_corner = ppy(ikmax);

% Locate the point on the discrete L-curve which is closest to the
% corner of the spline curve.  Prefer a point below and to the
% left of the corner.  If the curvature is negative everywhere,
% then define the leftmost point of the L-curve as the corner.
if (kmax < 0)
    reg_c = reg_param(lr); rho_c = rho(lr); eta_c = eta(lr);
else
    index = find(lrho < x_corner & leta < y_corner);
    if ~isempty(index)
        [~,rpi] = min((lrho(index)-x_corner).^2 + (leta(index)-y_corner).^2);
        rpi = index(rpi);
    else
        [~,rpi] = min((lrho-x_corner).^2 + (leta-y_corner).^2);
    end
    reg_c = reg_param(rpi); rho_c = rho(rpi); eta_c = eta(rpi);
end


end


