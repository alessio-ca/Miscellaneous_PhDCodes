function MSD = MSD_Maxwell_fun(t,Omega,Gamma,tau)
check = 0;
%Analytic espression for the MSD of a single-relaxation Maxwell fluid.
%The Log-sum-exp trick is used to ensure that no overflow occurs.

while(check < 2)
    %Factor statement
    exp1 = 2*t*Omega;
    exp2 = t*(Gamma + Omega);
    fac1 = 2*Gamma*Omega;
    fac2 = 3*Gamma^2*tau*Omega;
    fac3 = tau*Omega^3;
    fac4 = 4*t*Gamma*tau*Omega;
    first_sum = [log(Gamma^3*tau) + exp1, log(2*fac1) + exp2, log(fac2) + exp1, log(fac4*Gamma^2) + exp2,...
        log(3*Gamma*tau*Omega^2) + exp1, log(2*t*Omega^3) + exp2, log(fac3) + exp1];
    second_sum = [log(Gamma^2) + exp1, log(fac1) + exp1, log(2*t*Gamma^2*Omega) + exp2, log(2*fac2) + exp2,...
        log(Omega^2) + exp1, log(2*fac3) + exp2, log(fac4*Omega^2) + exp2];
    third_sum = Gamma^2 - Gamma^3*tau - fac1 + fac2 + Omega^2 - 3*Gamma*tau*Omega^2 + fac3;
    
    %Constant for the log-sum-exp trick
    first_c = max(first_sum);
    second_c = max(second_sum);
    check = check + 2;
    
    %Convert to arbitrary precision if factors are still too large
    if first_c > 1e6
        t=vpa(t);
        Omega=vpa(Omega);
        Gamma=vpa(Gamma);
        tau=vpa(tau);
        check = check - 1;
    end
end

%Log-sum-exp trick for the MSD expression
x1 = first_c + log(sum(exp(first_sum - first_c)));
x2 = second_c + log(sum(exp(second_sum - second_c)));
x3 = third_sum;
x4 = exp(x1 - x2) - 1;
x5 = log(x3) - exp2;
x6 = log(abs(x4)) - exp2 + x2;

MSD = (exp(x5) + sign(x4)*exp(x6))/(tau*Omega*(Omega^2 - Gamma^2)^2);
end