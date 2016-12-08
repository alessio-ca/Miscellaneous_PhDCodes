% fps=1000;
% folder1='C:\Users\ac2014\Documents\MATLAB\Git-MATLAB\Data_Microrheology\';
% cell_folder_list={'Simulation_1000FPS_Kel_1e-6_D_2e-6'};
% filename=['track_particle_output',char(cell_folder_list(1)),'.csv'];
% traj=dlmread([folder1,'Linked_',filename],',');
% traj(:,4)=[];
% traj=circshift(traj,[0 1]);

function [Corrx,Corry,Times]=CoarseGrainACF(datax,datay,fps,P,M,S)

%Calculate the auto-correlation function (ACF) of a data series (2D) using a
%multiple-tau correlator, a method for the calculation of time correlaton
%functions "on the fly" (Ramirez et al., 2010)

%Parameters:
%datax and datay are the data series of which you want to calculate the ACF
%fps is the frame rate of the acquisition system 
%P,M and S are the coarse-graining parameters (see Ramirez et al. (2010))

%Typical parameter choice:
%P=16;
%M=2;
%S=5;

if mod(P,M)~=0
    error('Ratio P/M is not an integer!')
end

%Array initialization
Levelx = zeros(S,P);
Corrx = zeros(S,P);
Accumx = zeros(S,1);

Levely = zeros(S,P);
Corry = zeros(S,P);
Accumy = zeros(S,1);

Count = zeros(S,P);
CountAccum = zeros(S,1);
Times=zeros(S,P);

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
            start=1;
        else
            start=P/M + 1;
        end
        
        Times(i,start:end) = (1/fps)*M^(i-1)*(start-1:P-1);
    end



%%%%%% Main loop %%%%%%

for i=1:length(datax)
    %Update base level with data points
    UpdateLevel(1,datax(i),datay(i));
    UpdateCorrCount(1);
    UpdateAccum(1,datax(i),datay(i));
    
    if mod(i,10000)==0
        disp(i)
    end
end

for i=1:S
    %Calculate times
    CalculateTimes(i);
end

%Normalized autocorrelations
Corrx=(Corrx./Count);
Corrx=Corrx./var(datax);

Corry=(Corry./Count);
Corry=Corry./var(datay);

%Discard empty levels & condense levels into a single array
id=~isnan(Corrx);
Times=Times(id);
Corrx=Corrx(id);
Corry=Corry(id);

%Apply correct sorting
[Times,id]=sort(Times);
Corrx=Corrx(id);
Corry=Corry(id);

end