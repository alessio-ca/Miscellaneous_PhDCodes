% fps=1000;
% folder1='C:\Users\ac2014\Documents\MATLAB\Git-MATLAB\Data_Microrheology\';
% cell_folder_list={'Simulation_1000FPS_Kel_1e-6_D_2e-6'};
% filename=['track_particle_output',char(cell_folder_list(1)),'.csv'];
% traj=dlmread([folder1,'Linked_',filename],',');
% traj(:,4)=[];
% traj=circshift(traj,[0 1]);

function [Corr,G]=CoarseGrainMicroRheo(datax,datay,fps,P,M,S)

%Calculate the G modulus of a data series (2D) using a
%multiple-tau correlator (a method for the calculation of time correlaton
%functions "on the fly" (Ramirez et al., 2010)) Fourier transform of a 
%function f(t) using the scheme by Evans et al. (2010) (which assumes that 
%f(t) vanishes for t<0 and is sampled at a finite set of data points (tk,fk) 
%which need not be equally spaced. It assumes 2D data (x(t) and y(t) series)).

%Parameters:
%datax and datay are the data series 
%fps is the frame rate of the acquisition system 
%P,M and S are the coarse-graining parameters (see Ramirez et al. (2010))


%Typical parameter choice:
%P=16;
%M=2;
%S=20;
if M < 2
    error('M must be > 2')
end
if M > P
    error ('M must be smaller than P')
end
if mod(P,M)~=0
    error('Ratio P/M is not an integer!')
end

%Array initialization
Levelx = zeros(S,P);
Corrx = zeros(S,P);
CorrScaledx = zeros(S,P);
Accumx = zeros(S,1);

Levely = zeros(S,P);
Corry = zeros(S,P);
CorrScaledy = zeros(S,P);
Accumy = zeros(S,1);

Count = zeros(S,P);
CountAccum = zeros(S,1);
Times=zeros(S,P);
Omega=zeros(S,floor(P/2));

%Data preparation for ACF calculation
datax=datax-mean(datax);
datay=datay-mean(datay);


%Subroutines:

%1. UpdateLevel:
%   Enter a new data points in the array i
    function UpdateLevel(i,x,y)
        %i is the level to be updated
        %x is the element to insert
        
        Levelx(i,2:P) = Levelx(i,1:P-1);
        Levelx(i,1) = x;
        Levely(i,2:P) = Levely(i,1:P-1);
        Levely(i,1) = y;
        
    end

%2. UpdateCorrCount:
%   Updates the corr function of level i
    function UpdateCorrCount(i)
        %i is the level to be updated
        
        if i == 1
            start=1;
        else
            start=P/M + 1;
        end
        
        Corrx(i,start:end) = Corrx(i,start:end) + Levelx(i,1).*Levelx(i,start:end);
        Corry(i,start:end) = Corry(i,start:end) + Levely(i,1).*Levely(i,start:end);
        Count(i,start:end) = Count(i,start:end) + 1;
        
        
    end

%3. UpdateAccum:
%   Updates the accumulator of level i
    function UpdateAccum(i,x,y)
        %i is the level to be updated
        %x is the element to insert
        Accumx(i) = Accumx(i) + x;
        Accumy(i) = Accumy(i) + y;
        CountAccum(i) = CountAccum(i) + 1;
        
        if CountAccum(i) == M
            %Push coarse-grained data to the next level
            UpdateLevel(i+1,Accumx(i)/M,Accumy(i)/M);
            UpdateCorrCount(i+1);
            UpdateAccum(i+1,Accumx(i)/M,Accumy(i)/M);
            %Reset counters of current level
            Accumx(i)=0;
            Accumy(i)=0;
            CountAccum(i)=0;
        end
    end

%4. CalculateTimes:
%   Calculate the time values of level i according to the coarse-graining procedure
    function CalculateTimes(i)
        %i is the level of interest
        if i == 1
            startT=1;
        else
            startT=P/M + 1;
        end
        
        Times(i,startT:end) = (1/fps)*M^(i-1)*(startT-1:P-1);
        Omega(i,1:end) = (0:floor(P/2)-1)*2*pi*fps/(P * M^(i-1));
    end

%5. FouEvansMulti:
%   Calculate the Fourier transform of the ACF following Evans et al.
%   (2010).
    function G=FouEvansMulti(i)
        %i is the last level achieved
        %Build full concatenated correlation function & vectors
        startT=P/M + 1;
        startO=floor((P/M)/2)+1;
        FullACF=zeros(P + (P - startT + 1)*(i-1),2);
        FullTimes=zeros(P + (P - startT + 1)*(i-1),1);
        FullOmega=zeros(floor(P/2) + (floor(P/2) - startO  + 1)*(i-1),1);
        G=zeros(floor(P/2) + (floor(P/2) - startO  + 1)*(i-1),2); 
        
        FullACF(1:P,1)= CorrScaledx(1,1:P);
        FullACF(1:P,2)= CorrScaledy(1,1:P);
        FullTimes(1:P)= Times(1,1:P);
        FullOmega(1:floor(P/2))= Omega(i,1:end);
        
        if i > 1
            for k = 2:i
                FullACF(P + 1  + (k-2) * (P-startT+1) : P + (k-1)*(P-startT+1),1) = CorrScaledx(k,startT:end);
                FullACF(P + 1  + (k-2) * (P-startT+1) : P + (k-1)*(P-startT+1),2) = CorrScaledy(k,startT:end);
                FullTimes(P + 1  + (k-2) * (P-startT+1) : P + (k-1)*(P-startT+1)) = Times(k,startT:end);
                FullOmega(1 + floor(P/2) + (k-2) * (floor(P/2)-startO+1): floor(P/2) + (k-1) * (floor(P/2)-startO+1))= Omega(i-k+1,startO:end);  
            end
        end
        %Vectors need to start at t>0
        FullACF(1,:)=[];
        FullTimes(1)=[];
        FullOmega(1)=[];
        G(1,:)=[];
        
        %Evans_fft is composed of three single terms and a summatory (see Evans et
        %al., (2010))
        
        %Calculate the single terms
        Cw = repmat(1i*FullOmega,1,2);
        Cw = Cw + [(1-exp(-1i*FullOmega*FullTimes(1)))*(FullACF(1,1)-1)/FullTimes(1),(1-exp(-1i*FullOmega*FullTimes(1)))*(FullACF(2,1)-1)/FullTimes(1)];
        
        %Calculate the summatory
        for i=1:length(FullOmega)
            Cw(i,:)=Cw(i,:)-[sum(diff(exp(-1i*FullOmega(i)*FullTimes)).*diff(FullACF(:,1))./diff(FullTimes)),sum(diff(exp(-1i*FullOmega(i)*FullTimes)).*diff(FullACF(:,2))./diff(FullTimes))];
        end
        
        %Split into real and imaginary parts
        %N.B! From now on, we have frequencies, not pulsations!
        Cw = [FullOmega,real(Cw),imag(Cw)];
        Cw(:,2:3) = Cw(:,2:3)./repmat(Cw(:,1).^2,1,2);
        Cw(:,4:5) = -(Cw(:,4:5)./repmat(-Cw(:,1).^2,1,2) + 1./repmat(Cw(:,1),1,2));
        Cw(:,6:7) = Cw(:,2:3).^2 + Cw(:,4:5).^2;
        
        %Calculate G moduli
        G=[Cw(:,1),-Cw(:,4:5)./(repmat(Cw(:,1),1,2).*Cw(:,6:7)),-Cw(:,2:3)./(repmat(Cw(:,1),1,2).*Cw(:,6:7))];        
    end

        
        



%%%%%% Main loop %%%%%%

for i=1:S
    %Calculate times
    CalculateTimes(i);
end

for i=1:length(datax)
    %Update base level with data points
    UpdateLevel(1,datax(i),datay(i));
    UpdateCorrCount(1);
    UpdateAccum(1,datax(i),datay(i));
    
    %Every 10000 steps, output ACF and G
    if mod(i,10000)==0
        disp(i)
        CorrScaledx=(Corrx./Count);
        CorrScaledx=CorrScaledx./var(datax);
        CorrScaledy=(Corry./Count);
        CorrScaledy=CorrScaledy./var(datay);
        
        LevelActive=ceil(log(i/P)/log(M)) + 1;
        G=FouEvansMulti(LevelActive);
    end
end

%At the end of the cyclr, recalculate ACF and G for the complete series
%Normalized autocorrelations
CorrScaledx=(Corrx./Count);
CorrScaledx=CorrScaledx./var(datax);

CorrScaledy=(Corry./Count);
CorrScaledy=CorrScaledy./var(datay);

LevelActive=ceil(log(i/P)/log(M)) + 1;
G=FouEvansMulti(LevelActive);




%Discard empty levels & condense levels of ACF into a single array
id=~isnan(CorrScaledx);
Times=Times(id);
CorrScaledx=CorrScaledx(id);
CorrScaledy=CorrScaledy(id);

%Apply correct sorting
[Times,id]=sort(Times);
Corrx=Corrx(id);
Corry=Corry(id);
Corr=[Times,Corrx,Corry];


end