function [k_corner,rho,eta,index,info] = CRIE_corner_New(rho,eta,fig)
%L_CORNER Locate the "corner" of the L-curve.
%
% [reg_c,rho_c,eta_c] =
%        l_corner(rho,eta,reg_param)
%        l_corner(rho,eta,reg_param,U,s,b,method,M)
%        l_corner(rho,eta,reg_param,U,sm,b,method,M) ,  sm = [sigma,mu]
%
% Locates the "corner" of the L-curve in log-log scale.
%
% It is assumed that corresponding values of || A x - b ||, || L x ||,
% and the regularization parameter are stored in the arrays rho, eta,
% and reg_param, respectively (such as the output from routine l_curve).
%
% If nargin = 3, then no particular method is assumed, and if
% nargin = 2 then it is issumed that reg_param = 1:length(rho).
%
% If nargin >= 6, then the following methods are allowed:
%    method = 'Tikh'  : Tikhonov regularization
%    method = 'tsvd'  : truncated SVD or GSVD
%    method = 'dsvd'  : damped SVD or GSVD
%    method = 'mtsvd' : modified TSVD,
% and if no method is specified, 'Tikh' is default.  If the Spline Toolbox
% is not available, then only 'Tikh' and 'dsvd' can be used.
%
% An eighth argument M specifies an upper bound for eta, below which
% the corner should be found.
% Per Christian Hansen, DTU Compute, January 31, 2015.

% Ensure that rho and eta are column vectors.
rho = rho(:); eta = eta(:);
[~,I] = sort(rho);
rho=rho(I);
eta=eta(I);

order = 4;  % Order of fitting 2-D spline curve.

% Initialization.
if (length(rho) < order)
    error('Too few data points for L-curve analysis')
end


if length(rho)~=length(eta)
    error('Vectors rho and eta must have the same length')
end
if length(rho)<3
    error('Vectors rho and eta must have at least 3 elements')
end

if (nargin < 3) | isempty(fig)
    fig = 0;        % Default is no figure.
elseif fig < 0,
    fig = 0;
end

info = 0;
fin = isfinite(rho+eta);    % NaN or Inf will cause trouble.
nzr = rho.*eta~=0;          % A zero will cause trouble.
kept = find(fin & nzr);

if isempty(kept)
    error('Too many Inf/NaN/zeros found in data')
end
if length(kept) < length(rho)
    info = info + 1;
    warning('Corner:warnings', ...
        ['Bad data - Inf, NaN or zeros found in data\n' ...
        '         Continuing with the remaining data'])
end
rho = rho(kept);            % rho and eta with bad data removed.
eta = eta(kept);
if any(rho(1:end-1)>rho(2:end)) || any(eta(1:end-1)<eta(2:end))
    info = info + 10;
    warning('Corner:warnings', 'Lack of monotonicity')
end
% Prepare for adaptive algorithm.
nP = length(rho);           % Number of points.
P = log10([rho eta]);       % Coordinates of the loglog L-curve.
V = P(2:nP,:)-P(1:nP-1,:);  % The vectors defined by these coordinates.
v = sqrt(sum(V.^2,2));      % The length of the vectors.
W = V./repmat(v,1,2);       % Normalized vectors.
clist = [];                 % List of candidates.
p = min(5, nP-1);           % Number of vectors in pruned L-curve.
convex = 0;                 % Are the pruned L-curves convex?
% Sort the vectors according to the length, the longest first.
[Y,I] = sort(v);
I = flipud(I);
% Main loop -- use a series of pruned L-curves. The two functions
% 'Angles' and 'Global_Behavior' are used to locate corners of the
% pruned L-curves. Put all the corner candidates in the clist vector.
while p < (nP-1)*2
    elmts = sort(I(1:min(p, nP-1)));
    
    % First corner location algorithm
    candidate = Angles( W(elmts,:), elmts);
    if candidate>0,
        convex = 1;
    end
    if candidate & ~any(clist==candidate)
        clist = [clist;candidate];
    end
    
    % Second corner location algorithm
    candidate = Global_Behavior(P, W(elmts,:), elmts);
    if ~any(clist==candidate)
        clist = [clist; candidate];
    end
    
    p = p*2;
end
% Issue a warning and return if none of the pruned L-curves are convex.
if convex==0
    k_corner = [];
    info = info + 100;
    warning('Corner:warnings', 'Lack of convexity')
    return
end
% Put rightmost L-curve point in clist if not already there; this is
% used below to select the corner among the corner candidates.
if sum(clist==1) == 0
    clist = [1;clist];
end
% Sort the corner candidates in increasing order.
clist = sort(clist);
% Select the best corner among the corner candidates in clist.
% The philosophy is: select the corner as the rightmost corner candidate
% in the sorted list for which going to the next corner candidate yields
% a larger increase in solution (semi)norm than decrease in residual norm,
% provided that the L-curve is convex in the given point. If this is never
% the case, then select the leftmost corner candidate in clist.
vz = find(abs(diff(P(clist,2))) ...  % Points where the increase in solution
    >= abs(diff(P(clist,1)))); % (semi)norm is larger than or equal
% to the decrease in residual norm.
if length(vz)>1
    if(vz(1) == 1),  vz = vz(2:end);  end
elseif length(vz)==1
    if(vz(1) == 1),  vz = [];  end
end
if isempty(vz)
    % No large increase in solution (semi)norm is found and the
    % leftmost corner candidate in clist is selected.
    index = clist(end);
else
    % The corner is selected as described above.
    vects = [P(clist(2:end),1)-P(clist(1:end-1),1) ...
        P(clist(2:end),2)-P(clist(1:end-1),2)];
    vects = sparse(diag(1./sqrt(sum(vects.^2,2)))) * vects;
    delta = vects(1:end-1,1).*vects(2:end,2) ...
        - vects(2:end,1).*vects(1:end-1,2);
    vv = find(delta(vz-1)<=0);
    if isempty(vv)
        index = clist(vz(end));
    else
        index = clist(vz(vv(1)));
    end
end
% Corner according to original vectors without Inf, NaN, and zeros removed.
k_corner = kept(index);

if fig  % Show log-log L-curve and indicate the found corner.
    figure(fig); clf
    diffrho2 = (max(P(:,1))-min(P(:,1)))/2;
    diffeta2 = (max(P(:,2))-min(P(:,2)))/2;
    loglog(rho, eta, 'k--o'); hold on; axis square;
    % Mark the corner.
    loglog([min(rho)/100,rho(index)],[eta(index),eta(index)],':r',... 
           [rho(index),rho(index)],[min(eta)/100,eta(index)],':r') 
    % Scale axes to same number of decades.
    if abs(diffrho2)>abs(diffeta2),
        ax(1) = min(P(:,1)); ax(2) = max(P(:,1));
        mid = min(P(:,2)) + (max(P(:,2))-min(P(:,2)))/2;
        ax(3) = mid-diffrho2; ax(4) = mid+diffrho2;
    else
        ax(3) = min(P(:,2)); ax(4) = max(P(:,2));
        mid = min(P(:,1)) + (max(P(:,1))-min(P(:,1)))/2;
        ax(1) = mid-diffeta2; ax(2) = mid+diffeta2;
    end
    ax = 10.^ax; ax(1) = ax(1)/2; axis(ax);
	xlabel('residual norm || A x - b ||_2')
	ylabel('solution (semi)norm || L x ||_2');
    title(sprintf('Discrete L-curve, corner at %d', k_corner));
end

% =========================================================================
% First corner finding routine -- based on angles
function index = Angles( W, kv)
    % Wedge products
    delta = W(1:end-1,1).*W(2:end,2) - W(2:end,1).*W(1:end-1,2);
    [mm kk] = min(delta);
    if mm < 0		% Is it really a corner?
        index = kv(kk) + 1;
    else			% If there is no corner, return 0.
        index = 0;
    end

% =========================================================================
% Second corner finding routine -- based on global behavior of the L-curve
function index = Global_Behavior(P, vects, elmts)
    hwedge = abs(vects(:,2));  % Abs of wedge products between
    % normalized vectors and horizontal,
    % i.e., angle of vectors with horizontal.
    [An, In] = sort(hwedge);   % Sort angles in increasing order.
    % Locate vectors for describing horizontal and vertical part of L-curve.
    count = 1;
    ln = length(In);
    mn = In(1);
    mx = In(ln);
    while(mn>=mx)
        mx = max([mx In(ln-count)]);
        count = count + 1;
        mn = min([mn In(count)]);
    end
    if count > 1
        I = 0; J = 0;
        for i=1:count
            for j=ln:-1:ln-count+1
                if(In(i) < In(j))
                    I = In(i); J = In(j); break
                end
            end
            if I>0, break; end
        end
    else
        I = In(1); J = In(ln);
    end
    % Find intersection that describes the "origin".
    x3 = P(elmts(J)+1,1)+(P(elmts(I),2)-P(elmts(J)+1,2))/(P(elmts(J)+1,2) ...
        -P(elmts(J),2))*(P(elmts(J)+1,1)-P(elmts(J),1));
    origin = [x3 P(elmts(I),2)];
    % Find distances from the original L-curve to the "origin".  The corner
    % is the point with the smallest Euclidian distance to the "origin".
    dists = (origin(1)-P(:,1)).^2+(origin(2)-P(:,2)).^2;
    [Y,index] = min(dists);


