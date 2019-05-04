function A = CRIE_Kernel(t,g,c,Kernel,Igrid)
%CRIE_Kernel Computes the matrix A from a Kernel and quadrature
%specification.
%
% CRIE_Kernel(t,g,c,Kernel)
%
% Given a time series t (m x 1), grid and quadrature vectors g and c (n x 1)
% and a Kernel specification, CRIE_Kernel computes the matrix A (m x n) for
% the problem Ax = b, with b a vector (m x 1) and x a vector (n x 1)
%
%
% FINISH DOCUMENTATION! Missing Igrid and Kernel description

% CREATED: Alessio Caciagli, University of Cambridge, April 2019

t = t(:); % must be column
g = g(:); % must be column
c = c(:); % must be column

if Igrid == 2
    g = log10(g);
end

%Kernel generation
[gM,~] = meshgrid(g,t);
[cM,tM] = meshgrid(c,t);

if Kernel == 1
    if Igrid == 2
        A = cM.*(log(10)*10.^(gM)).*exp(-tM.*10.^(gM));
    else
        A = cM.*exp(-tM.*gM);
    end
elseif Kernel == 2
    if Igrid == 2
        A = cM.*(log(10)*10.^(gM)).*(1 - exp(-tM.*10.^(gM)));
    else
        A = cM.*(1 - exp(-tM.*gM));
    end
else
    if all(size(Kernel) == [length(y),length(g)]) == 1
        if Igrid == 2
            A = cM.*(log(10)*10.^(gM)).*Kernel;
        else
            A = cM.*Kernel;
        end
    else
        error('Dimension mismatch. Check the dimensions of your Kernel.');
    end
end