function [s,g,yfit,lambda,info] = CRIE_New(t,y,Ng,Nl,input)
% CRIE Constrained Regularized method for inverting Intergal Equation (CRIE) 
% 
% [s,yfit,lambda,info] = CRIE(t,y,Ng,Nl,input) solves the integral equation
% 
%   y_k = int_a^b dg F_k(g) s(g) + sum_{i=1}^Nl L_{ki} b_i 
%        + sum_{i=1}^Nl L_{ki} beta_i 
%  
% where y_k denotes the observed data at k=1...Ny, F_k(g) represents the
% operator kernel, s(g) is the function to be estimated and the extra  
% optional sum (with beta_i unknowns) allows for background terms.
%
% The integral equation is discretized and solved according to the
% strategies and methods adopted in the CONTIN package. For information,
% see: http://www.s-provencher.com/contin.shtml.
% 
% The solution is obtained via a Regularized Least Squares (RLS)
% strategy. The regularization parameter is chosen accordingly to the
% L-curve method. 
% CRIE permits general linear equality and inequality constraints to be
% imposed on solution vector as: 
%   sum_{j+1}^Nx D_{ij} x_j >= d_i  i=1...Nnq
%   sum_{j+1}^Nx E_{ij} x_j == e_i  i=1...Neq
% where the matrices D and E, the arrays d and e and the variables Nnq and
% Neq can be specified.
%
% y: data array at discrete points t (for DLS analysis, y is the g1
% function and t represent the time points)
% Ng: number of grid points for the discretization and solution of the
% integral equation.
% Nl: number of term in the optional extra sum
%
% input is a structure containing the parameters for the calculation. The
% following parameters must be given:
%
%   - Iquad: specify the quadrature algorithm to be used 
%            trapezoidal = 1, simpson = 2
%   - Igrid: specify the grid in g-space to be used
%            linear = 1, logarithmic = 2
%   - Kernel: specify the kernel of the integral operator.
%             It's a Ny*Ng matrix.
%             If Kernel = 1, a Inverse Laplace transform kernel will be
%             used. If Kernel = 2, a Voigt fluid kernel will be used.
%   - Nnq: specify the number of inequalities
%   - Anq: specify the inequality matrix D and vector d.
%          It's a Nnq*(Nx+1) matrix. The first Nnq*Nx matrix specifies the
%          actual matrix D. The last column vector specifies d.
%   - Neq: specify the number of equalities
%          It's a Neq*(Nx+1) matrix. The first Neq*Nx matrix specifies the
%          actual matrix E. The last column vector specifies e.
%   - Nneg: specify the non-negativity condition  on s(g) (1=active, 0=inactive).
%           If active, it overrides Nnq and Anq.
%   - Ny0: specify the t=0 condition of y (1=active, 0=inactive).
%          If active, it overrides Neq and Aeq.
%   - iwt: specifies the weighted analysis of y.
%          unweighted = 1 , Poisson = 2, Squared = 3, UserInput = 4 
%
%   The following OPTIONAL parameters can be specified:
%   - g_lims: specify the upper and lower limits in the g-space. 
%             It's a 1x2 vector. Default: g_lims=[1/(3*t(end)),3/t(1)]
%   - wt: specifies the weights on y. It must be provided if iwt = 4.
%   - Nalpha: specifies the number of regularizer values probed.
%             Default = 40
%   - alpha_lims: specifies the upper and lower values of the probed
%                 regularizer values. Default: alpha_lims = [0.01,100]
%   - Nord: specifies the regularization order. 
%           1 = Ridge Regression, 2 = Second derivative (default)
%   - Nend: specifies the number of null points at the grid extrema. 
%           It's a 1x2 vector. Default: Nend = [2,2]
%   - KernelL: specifies the functions L_{ki} in the summation term.
%              It's a Ny*Nl matrix. 
%   - Nbg: specify the background condition (2=linear, 1=constant, 0=inactive).
%          If active, all functions L_{ki} in the summation term are constants.
%          It overrides KernelL.
%  
%
% s: the solution s(g) at the grid points g.
% g: the Ng grid points at which s(g) is calculated.
% yfit: fit of y_k according to the solution s and grid points g.
% lambda: optimal regularization parameter obtained via the L-curve.
%
% info is a structure containing the various parameters of the calculation.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% EXAMPLE: Inverse Laplace Transform of a noisy double-exponential decay
% function.
% % Data series:
% t = logspace(-7,0,200);
% y = 0.5*exp(-t/1e-3) + 0.5*exp(-t/1e-5);
% y = y + 0.07*(rand(size(t)) - 0.5);
% y(y>1) = 1;
% 
% % Parameters specification:
% Ng = 161;
% Nl = 0;
% ILT_Input.Iquad = 2;      % Simpson's rule
% ILT_Input.Igrid = 2;      % Log grid
% ILT_Input.Kernel = 1;     % ILT Kernel
% ILT_Input.Nnq = 0;        
% ILT_Input.Anq = 0;
% ILT_Input.Neq = 0;
% ILT_Input.Aeq = 0;
% ILT_Input.Nneg = 1;       % Non-negativity constraint
% ILT_Input.Ny0 = 1;        % y(0)=1 condition
% ILT_Input.iwt = 1;        % Unweighted analysis
% 
% [s,g,yfit,lambda,info] = CRIE(t,y,Ng,Nl,ILT_Input);

% CREATED: Alessio Caciagli, University of Cambridge, April 2018
%% Input
Ny = length(y);
if bitget(Ng,1) == 0 %(even)
    Ng = Ng + 1;
end
Nx = Ng + Nl;

Iquad = input.Iquad;
Igrid = input.Igrid;

%Grid points (if no grid is needed, set g_min = 1 and g_max = Ng)
%Include a fex extra grid points for numerical accuracy
%Default: g_min = 1/(3*t(end)), g_max = 3/t(1)
if isfield(input, 'g_lims')
    g_min = input.g_lims(1);
    g_max = input.g_lims(2);
else
    g_min = 1/(3*t(end));
    g_max = 1/(3*t(1));
end

%Kernel can be 1 (ILT), 2 (VF) or a Ny*Ng matrix
Kernel = input.Kernel;

% Specify inequalities/equalities.
% The matrix Anq/Aeq must be Nnq/Neq*(Ng+Nl+1) in size
Nnq = input.Nnq; %Number of inequalities
Anq = input.Anq; %Inequalities matrix
Neq = input.Neq; %Number of equalities
Aeq = input.Aeq; %Equalities matrix

% Specify non-negativity of solution
Nneg = input.Nneg;

% Specify first point constraint (y(0) = 1)
Ny0 = input.Ny0;

%Specify LS weights
iwt = input.iwt;
if isfield(input, 'wt')
    wt = input.wt;
else
    wt = ones(Ny,1);
end

% Specify regularizer settings
if isfield(input, 'Nalpha')
    Nalpha = input.Nalpha;
else
    Nalpha = 40;
end
if isfield(input, 'alpha_lims')
    alpha_min = input.alpha_lims(1);
    alpha_max = input.alpha_lims(2);
else
    alpha_min = 0.01;
    alpha_max = 100;
end

alphaV = logspace(log10(alpha_min),log10(alpha_max),Nalpha);

if isfield(input, 'Nord')
    Nord = input.Nord; %Order of regularizor
else
    Nord = 2;
end
if isfield(input, 'Nend')
    Nend = input.Nend; %Number of null external grid points
else
    Nend = [2,2];
end
if isfield(input, 'Nbg')
    Nbg = input.Nbg;
else
    Nbg = 0;
end

%% Grid generation
if Igrid == 1
    g = linspace(g_min,g_max,Ng);
elseif Igrid == 2
    g_min = log10(g_min);
    g_max = log10(g_max);
    g = linspace(g_min,g_max,Ng);
else
    error('Unknown grid specification. Check the parameter "Igrid".');
end
gm = g(1) - (g(2)-g(1));
gp = g(end) + (g(2)-g(1));
gpp = gp + (g(2)-g(1));

if Iquad == 1
    c = ones(1,Ng);
elseif Iquad == 2
    c = [4*ones(1,(Ng-3)/2);2*ones(1,(Ng-3)/2)]; %Simpson rule for coeff
    c = c(:);
    h = (g(end) - g(1))/(Ng-1);
    c = (h/3)*[1;c;4;2];
else
    error('Unknown quadrature specification. Check the parameter "Iquad".');
end

%% Input post-process
t = t(:); % must be column
y = y(:); % must be column
g = g(:); % must be column
c = c(:); % must be column

%Kernel generation
[gM,~] = meshgrid(g,t);
[cM,tM] = meshgrid(c,t);

if Kernel == 1
    if Igrid == 2
        A = cM.*(log(10)*10.^(gM)).*exp(-tM.*10.^(gM));
        Coeff = c.*(log(10)*10.^(g));
    else
        A = cM.*exp(-tM.*gM);
        Coeff = c;
    end
elseif Kernel == 2
    if Igrid == 2
        A = cM.*(log(10)*10.^(gM)).*(1 - exp(-tM.*10.^(gM)));
        Coeff = c.*(log(10)*10.^(g));
    else
        A = cM.*(1 - exp(-tM.*gM));
        Coeff = c;
    end
else
    if all(size(Kernel) == [Ny,Ng]) == 1
        if Igrid == 2
            A = cM.*(log(10)*10.^(gM)).*Kernel;
            Coeff = c.*(log(10)*10.^(g));
        else
            A = cM.*Kernel;
            Coeff = c;
        end
    else
        error('Dimension mismatch. Check the dimensions of your Kernel.');
    end
end

%Internal variables & equality/inequality constrains
UserL = 0;
if Nl~=0
    UserL  = 1;
end
UserNQ = 0;
if Nnq~=0
    UserNQ = 1;
    if all(size(Anq) == [Nnq,Nx+1]) ~= 1
        error('Dimension mismatch. Check the dimensions of your Anq.');
    end
end
UserEQ = 0;
if Neq~=0
    UserEQ = 1;
    if all(size(Aeq) == [Neq,Nx+1]) ~= 1
        error('Dimension mismatch. Check the dimensions of your Aeq.');
    end
end
if (Nneg ~= 0 && Nneg ~= 1)
    error('Non-negativity condition is not defined. Check Nneg.');
end

%Error estimation
%The first run is always unweighted unless the weights are user specified
if iwt == 4
    if length(wt) ~= Ny
        error('Error in user''s specification of the weights. Check wg');
    end
end

%Regularizor
Nreg = Ng - 2 + sum(Nend);
if Nord == 0
    R = eye(Ng,Ng);
    R(1,:) = [];
    R(end,:) = [];
    RUp = eye(Nend(1),max(Ng+Nend(1)-1,Ng));
    RUp(:,1:(Nend(1)-1)) = [];
    RDo = rot90(eye(Nend(2),max(Ng+Nend(1)-1,Ng)),2);
    RDo(:,(2 + end - Nend(2)):end) = [];
    R = [RUp;R;RDo];
elseif Nord == 2
    R = diag(-2 * ones(max(Nreg,Ng),1), 0) + diag(ones(max(Nreg,Ng)-1,1), -1) + diag(ones(max(Nreg,Ng)-1,1), 1);
    R(1:1-Nend(1),:) = [];
    R(end:end-Nend(2),:) = [];
    R(:,1:(Nend(1)-1)) = [];
    R(:,(2 + end - Nend(2)):end) = [];
else
    error('Unknown regularizor order specification. Check the parameter "Nord".');
end
R = [R,zeros(Nreg,Nl)]; %Do not regularize the constant term

%% Treatment of user input
if  UserL  == 1
    if Nbg == 2 %Linear background condition
        A = [A,repmat(t,1,Nl)]; %Add Nl columns of t's to the Kernel.
    elseif Nbg == 1 %Constant background condition
        A = [A,repmat(ones(Ny,1),1,Nl)]; %Add Nl columns of 1's to the Kernel.
    elseif Nbg == 0
        KernelL = input.KernelL; %Optional functions L_{ki}
        if all(size(KernelL) == [Ny,Nl]) == 1
            A = [A,KernelL];
        else
            error('Dimension mismatch. Check the dimensions of your KernelL.');
        end
    else
        error('Unknown Nbg specification. Check the parameter "Nbg".');
    end
else
    Nbg = 0;
end

%% Treatment of Equalities/Inequalities

% Equality constraints treatment (if needed)
if Ny0 == 1
    Neq = 1;
    if Igrid == 1
        Aeq = c';
    else
        Aeq = (c.*(log(10)*10.^(g)))'; 
    end
    if Nbg == 2
        Aeq = [Aeq,1,zeros(1,Nl)]; %y(0) = 1 constraint + bg terms (linear)
    elseif Nbg == 1
        Aeq = [Aeq,1,ones(1,Nl)]; %y(0) = 1 constraint + bg terms (constant)
    else
        Aeq = [Aeq,1]; %y(0) = 1 constraint
    end
end
Nxe = Nx - Neq;
% Inequality constaints treatment (like non-negativity)
% N.B! It holds also for the bg term!
if Nneg == 1
    Nnq = Nx;
    Anq = [eye(Nnq,Nx),zeros(Nnq,1)];
end
if Nreg < Nxe
    R = [R;zeros(Nxe-Nreg,Nx)];
end

%% Computation
Lx = zeros(size(alphaV));
sol = zeros(size(alphaV));
lambdaV = zeros(size(alphaV));

Aw = diag(sqrt(wt))*A;
yw = y.*sqrt(wt);

% Array scaling
%AsscaleM = mean(Aw,1);
%AsscaleS = std(Aw,1);
%Aw = (Aw - repmat(AsscaleM,size(Aw,1),1))./repmat(AsscaleS,size(Aw,1),1);


%QR factorization (if needed)
if Ny > Nx
    [C,eta]=CRIE_qr(Aw,yw);
else
    C = Aw;
    eta = yw;
end


options = optimoptions(@lsqlin,'Algorithm','active-set','Display','off','MaxIter',1500);
id ='optimlib:lsqlin:WillBeRemoved';
warning('off',id)

progressbar %Create progress bar

    
% Active-set solver with constraints
for i = 1:Nalpha
    alphaP=alphaV(i);
    if Aeq == 0
        if Anq == 0
            [s,~,~,exitflag] = lsqlin([C;alphaP*R],[eta(:);zeros(size(R,1),1)],[],[],[],[],[],[],[],options);
        else
            [s,~,~,exitflag] = lsqlin([C;alphaP*R],[eta(:);zeros(size(R,1),1)],-Anq(:,1:end-1),-Anq(:,end),[],[],[],[],[],options);
        end
    else
        [s,~,~,exitflag] = lsqlin([C;alphaP*R],[eta(:);zeros(size(R,1),1)],-Anq(:,1:end-1),-Anq(:,end),Aeq(:,1:end-1),Aeq(:,end),[],[],[],options);
    end
    sol(i) = norm(C*s - eta);
    Lx(i) = norm(R*s);
    lambdaV(i) = alphaP;
    if exitflag == -2
        sol(i) = [];
        Lx(i) = [];
        lambdaV(i) = [];
    end
    progressbar(i/Nalpha);
end
sol = nonzeros(sol);
Lx = nonzeros(Lx);
lambdaV = nonzeros(lambdaV);
% Corner of L_curve and recalculation with optimal alpha
if isempty(Lx)
    warning('CRIE:warnings', 'Optimal lambda has not been individuated. Selecting highest lambda in list.');
    index = length(lambdaV);
    k_c = index;
else
    [k_c,sol,Lx,index] =  CRIE_corner_New(sol,Lx);
end
if isempty(k_c)
    warning('CRIE:warnings', 'Optimal lambda has not been individuated. Selecting highest lambda in list.');
    index = length(lambdaV);
    k_c = index;
end
lambda = lambdaV(k_c);
if Aeq == 0
    if Anq == 0
        [s,~,~,~] = lsqlin([C;lambda*R],[eta(:);zeros(size(R,1),1)],[],[],[],[],[],[],[],options);
    else
        [s,~,~,~] = lsqlin([C;lambda*R],[eta(:);zeros(size(R,1),1)],-Anq(:,1:end-1),-Anq(:,end),[],[],[],[],[],options);
    end
else
    [s,~,~,~] = lsqlin([C;lambda*R],[eta(:);zeros(size(R,1),1)],-Anq(:,1:end-1),-Anq(:,end),Aeq(:,1:end-1),Aeq(:,end),[],[],[],options);
end

yfit = A*s;

%% Write to info structure
info.t = t;
info.y = y;
info.yfit = yfit;
info.grid = Igrid;
if Igrid == 1  
    info.g = g;
else
    info.g = 10.^g;
end
info.lambda = lambda;
info.c = c;
info.A = A;
info.Coeff = Coeff.*s(1:end-Nl);
info.Aeq = Aeq;
info.Anq = Anq;
info.R = R;
info.wt = wt;
info.date = datestr(now,30);
%% Plotting
close all

subplot(2,2,[1,3])
if Igrid == 1
    plot(t,y,'*',t,yfit)
else
    semilogx(t,y,'*',t,yfit)
end
xlabel('t')
ylabel('y(t)')

subplot(2,2,2)
if Igrid == 1
    plot(g,s(1:end - Nl,:),'-*')
else
    g = 10.^g;
    semilogx(g,s(1:end - Nl,:),'-*')
end
xlabel('s')
ylabel('g(s)')

subplot(2,2,4)
if ~isempty(Lx)
    loglog(sol,Lx,'k--o')
    hold on
    loglog([min(sol) - 0.2*(max(sol) - min(sol)),sol(index)],[Lx(index),Lx(index)],':r',...
        [sol(index),sol(index)],[0.5*min(Lx),Lx(index)],':r')
    hold off
    xlabel('residual norm || A x - b ||')
    ylabel('solution (semi)norm || L x ||');
end
end
