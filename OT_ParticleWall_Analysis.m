clear all
close all

HeightName = 'z=2d';

FileList = dir(['L:\OT_ParticleWall\*',HeightName,'*0.0025*.csv']);
minsize = 1e5;
Vx = zeros(minsize,2*length(FileList));
for i=1:length(FileList)
    filename = [FileList(i).folder,'\',FileList(i).name];
    data = dlmread(filename);
    if length(data) < minsize
        minsize = length(data);
        Vx(minsize+1:end,:)=[];
        t = data(:,1);
        Vx(:,(2*i - 1):2*i) = data(:,2:3);
    else 
        Vx(:,(2*i - 1):2*i) = data(1:minsize,2:3);
    end
end

%dlmwrite(['L:\OT_ParticleWall\',HeightName,'_G=0_0025.csv'],Vx,',')

Vx = Vx/10; %px to um conversion
fps = 1127;
[msd_w,tau_w,~,~] = msd_routine(1/fps,Vx);
[psd_w,f_w,~,~] = psd_routine(1/fps,Vx);

%%
FileList = dir(['L:\OT_ParticleWall\*',HeightName,'*0.0050*.csv']);
minsize = 1e5;
Vx = zeros(minsize,2*length(FileList));
for i=1:length(FileList)
    filename = [FileList(i).folder,'\',FileList(i).name];
    data = dlmread(filename);
    if length(data) < minsize
        minsize = length(data);
        Vx(minsize+1:end,:)=[];
        t = data(:,1);
        Vx(:,(2*i - 1):2*i) = data(:,2:3);
    else 
        Vx(:,(2*i - 1):2*i) = data(1:minsize,2:3);
    end
end

%dlmwrite(['L:\OT_ParticleWall\',HeightName,'_G=0_0050.csv'],Vx,',')

Vx = Vx/10; %px to um conversion
fps = 1127;
[msd_m,tau_m,~,~] = msd_routine(1/fps,Vx);
[psd_m,f_m,~,~] = psd_routine(1/fps,Vx);

%%
FileList = dir(['L:\OT_ParticleWall\*',HeightName,'*0.0100*.csv']);
minsize = 1e5;
Vx = zeros(minsize,2*length(FileList));
for i=1:length(FileList)
    filename = [FileList(i).folder,'\',FileList(i).name];
    data = dlmread(filename);
    if length(data) < minsize
        minsize = length(data);
        Vx(minsize+1:end,:)=[];
        t = data(:,1);
        Vx(:,(2*i - 1):2*i) = data(:,2:3);
    else 
        Vx(:,(2*i - 1):2*i) = data(1:minsize,2:3);
    end
end

%dlmwrite(['L:\OT_ParticleWall\',HeightName,'_G=0_0100.csv'],Vx,',')

Vx = Vx/10; %px to um conversion
fps = 1127;
[msd_s,tau_s,~,~] = msd_routine(1/fps,Vx);
[psd_s,f_s,~,~] = psd_routine(1/fps,Vx);

%%
figure(1)
loglog(tau_w,msd_w)
hold on
loglog(tau_m,msd_m)
loglog(tau_s,msd_s)
hold off
xlabel('Time [s]')
ylabel('MSD [\mum^2/s]')
legend('Weak','Medium','Strong')

figure(2)
loglog(f_w,psd_w)
hold on
loglog(f_m,psd_m)
loglog(f_s,psd_s)
hold off
xlabel('Frequency [Hz]')
ylabel('PSD [\mum^2/Hz]')
legend('Weak','Medium','Strong')

