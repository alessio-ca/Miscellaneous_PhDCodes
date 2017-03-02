function [MSD,Dxy,delay] = Langevin_MSD(M,dt,D)

% takes in the position vector M (1,n), the diffusion coefficient D and the time step dt
% returns a vector of MSDs for corresponding delay times MSD(t) = <x^2>(t) + <y^2>(t)
% returns the matrix Dxy(N,2), with Dx(t)=Dxy(t,1) and Dy(t)=Dxy(t,2),
% where Dx(t) and Dy(t) are the diffusion coefficients in x and y corresponding to delay time t

n=length(M(:,1));
times=(1:n)*dt;
MSD=zeros(n-1,1);
Dxy=zeros(n-1,2);
delay=zeros(n-1,1);
% if M has dimension n, then there will be n-1 different delay times
for j=1:n-1 % defines the step size btw data points in M
    delay(j)=j*dt; % save the corresponding delay time as the jth entry of delay
    sumr2=zeros(1,2);
    numbr2=0;
    k=1;
% in the following loop, we split up the trajectory into consecutive
% segments of length j in units of resolution time, starting at t=0.
% after each loop, we shift the whole chain of segments by one unit
% this is done as long as the total shift of the chain is smaller than the
% length of one segment and at least one segment can be fitted into the
% trajectory
    while k<j+1 && k+j<n+1
         Mnew=M(k:j:n,:);
         delta=diff(Mnew);
         r2=delta.^2;
         % compute the SD corresponding in x and y for each segment
         sumr2(1,1)=sumr2(1,1)+sum(r2(:,1));
         sumr2(1,2)=sumr2(1,2)+sum(r2(:,2));
         % add the sum of new SDs to the sum of existing SDs for step size j
         numbr2=numbr2+length(delta(:,1));
         % add the new number of SDs to the existing number of SDs
         k=k+1; % shift the chain by one unit
    end
    x2m=sumr2(1,1)/numbr2; % compute the overall MSD in x for step size j
    y2m=sumr2(1,2)/numbr2;% compute the overall MSD in y for step size j
    MSD(j)=x2m+y2m; % save this value as the jth entry of MSD
    Dxy(j,:)=[x2m/(2*delay(j)),y2m/(2*delay(j))]; % measured diffusion coefficients in x and y, for corresponding delay time
end
plot(delay,MSD,'.',times,4*D*times,'k-');
title('MSD vs Delay Time','FontSize',20);
xlim([0,0.2]);
xlabel('delay time (s)','FontSize',20);
ylabel('MSD (m2)','FontSize',20);
legend('Langevin sim.','4Ddt');
end