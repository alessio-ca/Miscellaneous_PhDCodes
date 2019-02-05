%This routine corrects a single-particle tracked trajectories for drift
%effects. 

clear all
close all

%Input: enter the last stamp of the track_output files (a list if more than
%one is present)
[cell_folder_list,fullpath]=uigetfile('*.csv','Select a file(s)','MultiSelect','on');
cell_folder_list=strrep(cell_folder_list,'.csv','');

if length(cell_folder_list)==1
    if cell_folder_list==0
        error('Canceled operation.')
    end
end
if ~iscell(cell_folder_list)
    cell_folder_list={cell_folder_list};
end
listlen=length(cell_folder_list);
%%
%Drift parameters:
% drift - order of polynomial correction (default is 1)
drift = 1;
%%
%Routine going
for frame=1:listlen
    filename=[char(cell_folder_list(frame)),'.csv'];
    disp(char(cell_folder_list(frame)))
    
    
    data=dlmread([fullpath,filename],',');
    positions=data(:,[2:3,1]);
    disp(['Substracting drift using ',num2str(drift),' order polynomial...']);
    
    %Substract mean
    positions(:,1:2) = positions(:,1:2) - mean(positions(:,1:2)); 
    datacopy = positions;
    
    %Moving average to "clean" trajectory
    Nd = 1000;
    if Nd ~= -1 && (Nd < 2 || Nd > length(positions))
            Nd = floor(length(linked)/500);
            disp(['Invalid Nd! New value is set to: ',num2str(Nd)]);
    end
    datacopy(:,1:2) = movmean(positions(:,1:2),Nd);
   
    %Fit n-order polynomial
    tfit = positions(:,3) - mean(positions(:,3));
    px=polyfit(tfit,datacopy(:,1),1);
    py=polyfit(tfit,datacopy(:,2),1);
    
    %Substract mean
    positions(:,1)=positions(:,1)-polyval(px,tfit);
    positions(:,2)=positions(:,2)-polyval(py,tfit);

    figure(1)
    subplot(1,2,1)
    plot(data(:,1),data(:,2) - mean(data(:,2)),'b')
    hold on
    plot(datacopy(:,3),datacopy(:,1),'k--') 
    plot(datacopy(:,3),polyval(px,tfit),'r.-')
    plot(positions(:,3),positions(:,1),'g')
    hold off
    title(['X Coordinate over time, series ',num2str(frame)])
    xlabel('Time step')
    ylabel('Pixel')
    legend('Raw','Rolling Mean','Polynomial fit','Dedrifted data')
    
    subplot(1,2,2)
    plot(data(:,1),data(:,3) - mean(data(:,3)),'b')
    hold on
    plot(datacopy(:,3),datacopy(:,2),'k--') 
    plot(datacopy(:,3),polyval(py,tfit),'r.-')
    plot(positions(:,3),positions(:,2),'g')
    hold off
    title(['Y Coordinate over time, series ',num2str(frame)])
    xlabel('Time step')
    ylabel('Pixel')
    legend('Raw','Rolling Mean','Polynomial fit','Dedrifted data')
    
    dlmwrite([fullpath,'De-drifted_',filename],positions(:,[3,1:2]), 'delimiter', ',', 'precision', 9);
    disp('Printed position trajectories.');

    
    figure(2)
    x=positions(:,1);
    y=positions(:,2);
    h=histogram(x,50,'Normalization','pdf');
    h_MidPoints = 0.5*(h.BinEdges(2:end) + h.BinEdges(1:end-1));
    hold on
    k=histogram(y,50,'Normalization','pdf');
    k_MidPoints = 0.5*(k.BinEdges(2:end) + k.BinEdges(1:end-1));
    hold off
    ylabel('Count')
    xlabel('Position (um)')
    title(['Histo data set ',num2str(frame,'%u')])
    drawnow
    fid = fopen([fullpath,'Histogram_',filename], 'w');
    fprintf(fid, ['x_Bins,' 'x_Values,' 'y_Bins,' 'y_Values\n']);
    fclose(fid);
    
    dlmwrite([fullpath,'Histogram_',filename],[h_MidPoints',h.Values',k_MidPoints',k.Values'],'-append', 'precision', 9);
    disp('Printed position histogram.');
    
    
    
    clear data
    clear datacopy
    clear positions
    fprintf('\n')
end

