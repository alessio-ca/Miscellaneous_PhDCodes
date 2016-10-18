%This routine executes a single-particle tracking for bright field images
%in which the particle has a doughnut-like appearance. It applies an
%2D autocorrelation filter with a custom filter input by the user and
%paraboloid fitting to get sub-pixel resolution (see AnalyseTrack_Parallel.m
%documentation for more details).
%It DOESN'T require pre-tracking procedures.
%It support parallelization.

%Input: enter mother directory name (it must end with a "\" sign)
folder1='/Users/AleCac/Desktop/Video_JAvi/18_Aug/';


%Input: enter child directory name
folder2='T=50C.18Aug2016_23.58.44';

%Input: enter filename stamp
stamp='T=50C.18Aug2016_23.58.44_';


img=imread([folder1,folder2,'/',stamp,'000000.tiff']); %Format this if you series doesn't start with 0 or has a different number of digits
imshow(mat2gray(img));
[filter{1},r_max]=CreateFilter(img);
close all

listlen = length(dir([folder1,folder2,'/','*.tiff']));
%listlen=804886;
%Input: by default, the program process all the pictures. Modify the 2nd
%and 3rd field if necessary
[poslist,histostep]=AnalyseTrack_Parallel([folder1,folder2,'/',stamp],1,listlen-1,filter,r_max*2); %The listlen-1 accounts for one frame being 00000
figure(1)
histogram(histostep);


%Input: by default, the program will create an output file in the mother
%directory. Modify if necessary
fullFileName=fullfile(folder1,['parallel_track_particle_output',folder2,'.csv']);
if exist(fullFileName,'file')==2
    %Input to specify how to deal with pre-existing file
    loop='A';
    while (loop ~= 'Y') && (loop ~= 'N')
        loop=input('Output track file already exists. Do you want to overwrite it? [Y/N] ','s');
        if loop=='Y'
            delete(fullFileName);
        elseif loop=='N'
            while exist(fullFileName,'file')==2
                fullFileName = [fullFileName(1:end-4),'bis.csv'];
            end
        else
            disp('Incorrect input from user.')
            loop = 'A';
        end
    end
end

dlmwrite(fullFileName, poslist, 'delimiter', ',', 'precision', 9); 
%%
figure(2)
plot(poslist(:,1),poslist(:,2))
hold on
plot(poslist(:,1),poslist(:,3))
hold off
title('Coordinates over time')
xlabel('Time step')
ylabel('Pixel')
legend('X coordinate','Y coordinate')

