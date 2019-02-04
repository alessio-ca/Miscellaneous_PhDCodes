P = @(q,R) (3*(R.^3).*(sin(q.*R) - q.*R.*cos(q.*R))./(q.*R).^3).^2;  %Single particle form factor

q = logspace(-2,0,1000)'; %q-space

%Simpson Integration
Ng = 1001; %Grid Integration
R = linspace(50,130,Ng)'; %R-space
c = [4*ones(1,(Ng-3)/2);2*ones(1,(Ng-3)/2)]; %Simpson rule for coeff
c = c(:);
h = (R(end) - R(1))/(Ng-1);
c = (h/3)*[1;c;4;1];
[cM,~] = meshgrid(c,q);
[RM,qM] = meshgrid(R,q);
Kernel = P(qM,RM);
A = cM.*Kernel;

Rm = 90;
Rstd1 = 2;
Rstd2 = 5;
Rd1 = (1/sqrt(2*pi*Rstd1))*exp(-(R - Rm).^2/(2*Rstd1^2)); %R distribution (mean=Rm,std=Rstd1)
Rd2 = (1/sqrt(2*pi*Rstd2))*exp(-(R - Rm).^2/(2*Rstd2^2)); %R distribution (mean=Rm,std=Rstd2)
Rd3 = (1/sqrt(2*pi*Rstd1))*exp(-(R - Rm).^2/(2*Rstd1^2))*0.8 + ...
       (1/sqrt(2*pi*Rstd1))*exp(-(R - Rm - 15).^2/(2*Rstd1^2))*0.2;

norm = 1/Rm^6;
y0 = norm*P(q,Rm); %P(q) for monomodal distro 
y1 = norm*A*Rd1; %P(q) from distro 1
y2 = norm*A*Rd2; %P(q) from distro 2
y3 = norm*A*Rd3; %P(q) from distro 2


loglog(q,y0,q,y1,q,y2,q,y3)
xlabel('q(nm)')
ylabel('P(q)')
legend('Monomodal','Distro1','Distro2')



