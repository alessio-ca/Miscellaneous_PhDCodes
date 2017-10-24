 n = 1.33;
tail = [1e-1,1];
%% Constants
kB = 1.38*10^(-23);
lambda = 685*10^(-9); %Laser wavelength
L = 2e-3; %Cuvette thickness
l_star = 401.87e-6; %Transport mean free path
k0 = 2*pi*n/lambda;
R = 115e-9; %Bead radius

%% g1(tau) calc. and multiexponential fit to g1(tau) (for smoothness)

alpha = 1e-2;
I_vec = zeros(100,size(data,2));

for j = 1:size(data,2)
    
    %Delete first 2 points (usually noise) & format data
    Time_ToFit = t(3:end);
    ACF_ToFit = data(3:end,j);
    if all(tail)
        tailFit = ACF_ToFit(Time_ToFit > tail(1) & Time_ToFit < tail(2));
        tailFit = mean(tailFit);
        ACF_ToFit = ACF_ToFit - tailFit;
    end
    CritPts = find(ACF_ToFit< 0.1*ACF_ToFit(1));
    tTemp = Time_ToFit(1:CritPts(1));
    
    %g1(tau) calculation with  optional tail treatment (least-Squares)
    if all(tail)
        dataTempLog = log(ACF_ToFit(1:CritPts));
        A = [ones(length(tTemp),1),tTemp,tTemp.^2];
        Coeff = A\dataTempLog;
        beta = exp(Coeff(1));
    else
        beta = 1;
    end
    dataTemp = sqrt(ACF_ToFit(1:CritPts))/sqrt(beta);
    
    U = tTemp;
    W = tTemp.^2;
    Z = dataTemp.*tTemp;
    V = dataTemp.*tTemp.^2;
    F = dataTemp.*tTemp.^3;
    H = dataTemp.*tTemp.^4;
    K = dataTemp;

    
    A = [W,U,ones(length(U),1),-V,-Z,-K];
    Coeff = A\F;
    strange_fit_23 = @(x) (Coeff(1).*x.^2 + Coeff(2).*x + Coeff(3))./(x.^3 + Coeff(4).*x.^2 + Coeff(5).*x + Coeff(6));
    fo = fitoptions('rat23');
    fo.StartPoint = Coeff';
    
    A = [U,ones(length(U),1),-V,-Z,-K];
    Coeff = A\F;
    strange_fit_13 = @(x) (Coeff(1).*x + Coeff(2))./(x.^3 + Coeff(3).*x.^2 + Coeff(4).*x + Coeff(5));
    go = fitoptions('rat13');
    go.StartPoint = Coeff';
    ro = fitoptions('rat33');
    ro.StartPoint = [0,0,Coeff'];
    
    A = [W,U,ones(length(U),1),-F,-V,-Z,-K];
    Coeff = A\H;
    strange_fit_24 = @(x) (Coeff(1).*x.^2 + Coeff(2).*x + Coeff(3))./(x.^4 + Coeff(4).*x.^3 + Coeff(5).*x.^2 + Coeff(6).*x + Coeff(7));
    ho = fitoptions('rat24');
    ho.StartPoint = Coeff';
    
    
    
    f = fit(tTemp, dataTemp, 'rat23',fo);
    fit1 = feval(f,Time_ToFit);
    g = fit(tTemp, dataTemp, 'rat13',go);
    fit2 = feval(g,Time_ToFit);
    h = fit(tTemp, dataTemp, 'rat24',ho);
    fit3 = feval(h,Time_ToFit);
    l = fit(tTemp, dataTemp, 'rat33',ro);
    fit4 = feval(l,Time_ToFit);

    
    semilogx(Time_ToFit,sqrt(ACF_ToFit/beta),'x')
    hold on
    semilogx(Time_ToFit,fit1,'b')
    semilogx(Time_ToFit,strange_fit_23(Time_ToFit),'b--')
    semilogx(Time_ToFit,fit2,'r')
    semilogx(Time_ToFit,strange_fit_13(Time_ToFit),'r--')
    semilogx(Time_ToFit,fit3,'g')
    semilogx(Time_ToFit,strange_fit_24(Time_ToFit),'g--')
    semilogx(Time_ToFit,fit4,'k')
    hold off
    

end