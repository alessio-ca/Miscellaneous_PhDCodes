%Input: enter the last stamp of the track_output files (a list if more than
%one is present)
fullpath='/Users/AleCac/Documents/MATLAB/Git-MATLAB/Data_Microrheology/';
cell_folder_list={'60x_1.20_Si02_1000FPS_1.5um_Tweezer_G_0.05.13Aug2016_15.27.31'};

listlen=length(cell_folder_list);

%Input: enter frame rate and pixel-to-nm conversion factor
fps=1000;
umPerPixel=0.1425;

%Input: drift correction order (-1 for no drift correction)
drift=1;
%Input: MSD calculation (1 for yes, 0 for no)
MSD=0;
%Input: Allan variance calculation (1 for yes, 0 for no)
allvar=1;

close all
clear pos
pos=cell(1,listlen);
%%
%Routine going
for frame=1:listlen
    %Param explanation (see documentation of file 'track.m'):
    %   'mem' = number of frames allowed for skipping
    %   'good' = minimum length of trajectory allowed
    %   'dim' = dimension (2D or 3D tracking)
    %   'quite' = 1 means screen output suppressed, 0 means screen output allowed
    
    %N.B! File containing particle positions must be formatted as
    %'track_particle_output[name_of_folder].csv'. If this is not the case,
    %modify the dlmread command below accordingly (even better: CHANGE YOUR
    %FILENAME AND STANDARDIZE YOURSELF, FOOL!)
    param=struct('mem',0,'good',100,'dim',2,'quiet',1);
    filename=['track_particle_output',char(cell_folder_list(frame)),'.csv']; 
    try linked=dlmread([fullpath,'Linked_',filename],',');
       disp('Linked trajectory file found.')
       
    catch ME
        if strcmp(ME.identifier,'MATLAB:dlmread:FileNotOpened')
            disp('No linked trajectory file found. Proceeding with linking...')
            data=dlmread([fullpath,filename],',');
            positions=data(:,[2:3,1]);            
            
            linked=track(positions,10,param);
            
            %Data Analysis is done only on the longest track
            maxidx=max(linked(:,4));
            maxidx=(1:maxidx);
            for i=1:length(maxidx)
                maxidx(i)=length(linked(linked(:,4)==maxidx(i)));
            end
            [~,maxidx]=max(maxidx);
            linked=linked(linked(:,4)==maxidx,:);
            
            linked(:,1:2)=linked(:,1:2)*umPerPixel;
            linked(:,3)=linked(:,3)*(1/fps);
            dlmwrite([fullpath,'Linked_',filename],linked, 'delimiter', ',', 'precision', 9);
        end
    end
    if drift==-1
        disp('No drift substraction. Substracting offset...');
        const=linked(1,1:2);
        linked(:,1:2)=linked(:,1:2)-repmat(const,length(linked),1);
    elseif drift >= 0
        disp(['Substracting drift using ',num2str(drift),' order polynomial...']);
        Nd=1000;
        if Nd ~= -1 && (Nd < 2 || Nd > length(linked))
            Nd = floor(length(linked)/500);
            disp(['Invalid Nd! New value is set to: ',num2str(Nd)]);
        end
        linkedcopy=linked;
        
        if Nd > 0
            for i=1:length(linked)-Nd
                linkedcopy(i,1:2) = mean(linked(i:i+Nd,1:2));
            end
            for i=0:Nd
                linkedcopy(end-i,1:2)= mean(linked(end-i-Nd:end-i,1:2));
            end
        end
        tfit = linked(:,3) - mean(linked(:,3));
        px=polyfit(tfit,linkedcopy(:,1),1);
        py=polyfit(tfit,linkedcopy(:,2),1);
        linked(:,1)=linked(:,1)-polyval(px,tfit);
        linked(:,2)=linked(:,2)-polyval(py,tfit);
        if ~exist([fullpath,'De-drifted_',filename],'file')
            dlmwrite([fullpath,'De-drifted_',filename],linked(:,[3,1:2]), 'delimiter', ',', 'precision', 9);
        end
    else
        error('Wrong drift parameter!')
    end
    
    figure(1)
    x=linked(:,1);
    y=linked(:,2);
    h=histogram(x,50,'Normalization','pdf');
    h_MidPoints = 0.5*(h.BinEdges(2:end) + h.BinEdges(1:end-1));
    hold on
    k=histogram(y,50,'Normalization','pdf');
    k_MidPoints = 0.5*(k.BinEdges(2:end) + k.BinEdges(1:end-1));
    hold off
    ylabel('Count')
    xlabel('Position (um)')
    title(['Histo data set ',num2str(frame,'%u')])
    fid = fopen([fullpath,'Histogram_',filename], 'w');
    fprintf(fid, ['# x_Bins\t' 'x_Values\t' 'y_Bins\t' 'y_Values\n']);
    fclose(fid);
    dlmwrite([fullpath,'Histogram_',filename],[h_MidPoints',h.Values',k_MidPoints',k.Values'],'-append', 'precision', 9);
    disp('Printed position histogram.');

    if MSD == 1
        disp('Calculating MSD...')
        pos{frame}=MSD_single_track(linked);
        pos{frame}(:,1)=pos{frame}(:,1)*1/fps;
        pos{frame}=pos{frame}(1:floor(length(pos{frame})/4),:);
        dlmwrite([fullpath,'MSD_',filename],pos{frame}, 'delimiter', ',', 'precision', 9);
    elseif MSD == 0
         try pos{frame}=dlmread(['MSD_',filename],',');
            disp('MSD file found.')
         catch ME
            if strcmp(ME.identifier,'MATLAB:dlmread:FileNotOpened')
                disp('MSD file not found! MSD will not be displayed.')
            else
                rethrow(ME)
            end
         end
    else
        error('Wrong MSD parameter!')
    end
    
    if allvar == 1
        disp('Calculating Allan variances...')
        avar{frame}=allan(linked(:,1),fps,logspace(log10(1/fps),log10((1/fps)*length(linked(:,1)))));
        dlmwrite([fullpath,'Allan_',filename],avar{frame}, 'delimiter', ',', 'precision', 9);
    elseif allvar == 0
        try avar{frame}=dlmread(['Allan_',filename],',');
            disp('Allan file found.')
        catch ME
            if strcmp(ME.identifier,'MATLAB:dlmread:FileNotOpened')
                disp('Allan file not found! Allan will not be displayed.')
            else
                rethrow(ME)
            end
        end
    else
        error('Wrong Allan parameter!')
    end
    
    
    
    clear data
    clear linkedcopy
    clear linked
end
%%
%Plot section
if MSD==1
    figure(2)
    hold on
    for frame=1:listlen
        x=pos{frame}(:,1);
        y=pos{frame}(:,2);
        n=pos{frame}(:,3);
        loglog(x,y);
        ylabel('MSD (um^{2})')
        xlabel('Time (s)')
    end
    hold off
    title(['MSD plot'])
    set(gca,'xscale','log');
    set(gca,'yscale','log');
end
%%
if allvar==1
    figure(3)
    hold on
    for frame=1:listlen
        x=avar{frame}(:,1);
        y=avar{frame}(:,2);
        loglog(x,y,'o');
        ylabel('Allan variance (um^{2})')
        xlabel('Time (s)')
    end
    hold off
    title(['Allan variance plot'])
    set(gca,'xscale','log');
    set(gca,'yscale','log');
end
