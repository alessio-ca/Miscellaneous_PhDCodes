function [x,point_out] = graddesc(x0, gradf_x, gradf_y, c0, r0, xmax, xmin, varargin)
%GRADDESC Gradient descent optimization.
%

% Point output
pt_out = 0;
for n = 1:2:length(varargin)
    if strcmpi(varargin{n},'PointOut')
        pt_out = varargin{n+1};
    end
end
if pt_out~=0 && pt_out ~=1
    error('Error: PointOut must be 0 or 1')
end

% Check function string
gradf_x = fcnchk(gradf_x);
gradf_y = fcnchk(gradf_y);



% Steepest descent parameters

x = x0;
niters = 200;
alpha0= 500;
dxmin = 1e-12; %Perturbation tolerance
tol = 1e-18; %Gradient tolerance;


if pt_out==1
    point_out = zeros(niters,2);
else
    point_out=0;
end

chkcircles = zeros(size(c0,1));
%  Main optimization loop.
alpha = alpha0;
for j = 1:niters
  grad_x = gradf_x(x');
  grad_y = gradf_y(x');
  dx = alpha*[grad_x,grad_y];
  x =  x + dx;
  if pt_out==1
      point_out(j,:)=x;
  end
  
  gnorm = norm([grad_x,grad_y]);
  dxnorm = norm(dx);
  chkcircles = sum((x-c0).^2,2) < r0^2;
  
  if (gnorm < tol || dxnorm < dxmin || ~all(xmax-x>0) || ~all(x-xmin>0) || any(chkcircles))
    % Termination criteria are met
    point_out(j+1:end,:)=repmat(x,[length(point_out)-j,1]);
    if  any(chkcircles)
        disp('Particle is in contact with another.')
    elseif ~all(xmax-x>0) || ~all(x-xmin>0) 
        disp('Particle is outside the ROI.')
    else
        disp('Equilibrium position reached')
    end
    
    x(3)=1;
    return;
  end
  dgrad = -[gradf_x(x'),gradf_y(x')] + [grad_x,grad_y];
  alpha = (dx * dgrad')./norm(dgrad)^2;
  
end

if (j == niters)
  disp('Warning: Maximum number of iterations has been exceeded in graddesc');
  x(3)=0;
end