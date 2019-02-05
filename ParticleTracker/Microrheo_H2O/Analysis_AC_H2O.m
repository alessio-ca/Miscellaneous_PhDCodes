close all
clear all

FileList = dir('T*');
Legend = cell(size(FileList));
subplot(1,2,1)
hold on
for i=1:length(FileList)
    filename = [FileList(i).folder,'\',FileList(i).name];
    load(filename)
    Legend{i} = FileList(i).name(3:4);
    if i <= 6
        plot(tau,ACF,'--o')
    else
        plot(tau,ACF,'--s')
    end
end
hold off
set(gca,'XScale','log')
legend(Legend)
xlabel('Time (us)')
ylabel('g1')


n=1.33;
theta = 173*2*pi/360;
lambda = 633*10^(-9);
q = 4*pi*n*sin(theta/2)/lambda;
subplot(1,2,2)
hold on
for i=1:length(FileList)
    filename = [FileList(i).folder,'\',FileList(i).name];
    load(filename)
    MSD = 1e+6*(6/q^2)*(-log(ACF));
    Legend{i} = FileList(i).name(3:4);
    if i <= 6
        plot(tau,MSD,'--o')
    else
        plot(tau,MSD,'--s')
    end
end
hold off
set(gca,'XScale','log')
set(gca,'YScale','log')
legend(Legend,'Location','southeast')
xlabel('Time (us)')
ylabel('MSD (m^2/s)')

