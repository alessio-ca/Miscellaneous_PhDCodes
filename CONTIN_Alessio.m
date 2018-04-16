% CREATED: Alessio Caciagli, University of Cambridge, March 2018
clear variables;close all

%Fictitious data
t = logspace(-7,2,200);
y = exp(-t/1e-3);

%% Input
Ny = length(y); 
Ng = 161; %Number of grid points
Nl = 0; %Number of optional terms
Nx = Ng + Nl;
if bitget(Ng,1) == 0 %(even)
    Ng = Ng + 1;
end

Iquad = 1; %direct solution
Iquad = 2; %quadrature solution

Igrid = 1; %linear grid
Igrid = 2; %log grid
%Grid points (if no grid is needed, set s_min = 1 and s_max = Ng)
%N.B! Remember to include a fex extra grid points for numerical accuracy!!!
g_min = 1/(3*t(end));
g_max = 3/t(1);

%Kernel can be 1 (ILT) or a Ny*Ng matrix
Kernel = 1; %Inverse Laplace Transform kernel

% Specify inequalities/equalities.
% The matrix Anq/Aeq must be Nnq*(Ng+Nl+1) in size 
Nnq = 0; %Number of inequalities
Anq = 0;
Neq = 2; %Number of equalities
Aeq = randi(10,Neq,Nx+1);

% Specify non-negativity of solution
Nneg = 0; %No constains
Nneg = 1; %Constrain to be positive

%Specify LS weights
iwt = 1; %Unweighted analysis
iwt = 2; %Poisson's stat 
iwt = 3; %Constant relative error
iwt = 4; %User specified (in wt)
wt = ones(Ny,1);

% Specify regularizer settings
Nalpha = 20;
alphaV = logspace(-4,2,Nalpha);
Nord = 2; %Order of regularizor
Nend = [2,2]; %Number of null external grid points

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
%g0 = g0(:); % must be column

%Kernel generation
[gM,~] = meshgrid(g,t);
[cM,tM] = meshgrid(c,t);
if Kernel == 1
    if Igrid == 2
        A = cM.*(log(10)*10.^(gM)).*exp(-tM.*10.^(gM));
    else
        A = cM.*exp(-tM.*gM);
    end
else
    if all(size(Kernel) == [Ny,Ng]) == 1
        A = Kernel;
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
switch iwt
    case 4
        if length(wt) ~= Ny
            error('Error in user''s specification of the weights. Check wg');
        end
    otherwise
        wt = ones(Ny,1);
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
%R = [[1,zeros(1,Ng-1)];diag(-2 * ones(Ng,1), 0) + diag(ones(Ng-1,1), -1) + diag(ones(Ng-1,1), 1);[zeros(1,Ng-1),1]];
%R = (1/h)^2 * R;
%% Computation
s=zeros(size(g));
% Array rescaling to interval [0,1]

%QR factorization (if needed)
if Ny > Nx
    [A,y]=CONTIN_Alessio_qr(diag(1./sqrt(wt))*A,y./sqrt(wt));
end

% Equality constraints treatment (if needed)
if Neq > 0
    K = CONTIN_Alessio_Q(Aeq(:,1:end-1)');
    K1 = K(:,1:Neq);
    K2 = K(:,Neq+1:end);
    s1 = (Aeq(:,1:end-1)*K1)\Aeq(:,end);
else
    K1 = [];
    s1 = [];
    K2 = eye(Nx);
end
Nxe = Nx - Neq;
if size(R,1) < Nxe
    R = [R;zeros(Nxe - size(R,1),size(R,2))];
end

% Inequality constaints treatment (like non-negativity)
if Nneg == 1
    Nnq = Ng;
    Anq = [eye(Nnq,Nx),zeros(Nnq,1)];
end
if Nreg < Nxe
    R = [R;zeros(Nxe-Nreg,Ng)];
end

%SVD decomposition
[U,H1,Z] = svd(R*K2);
H2 = H1(1 + end - (Nreg - Nxe):end,:);
H1 = H1(1:Nxe,:);
    
    


    
    

        
        
        




