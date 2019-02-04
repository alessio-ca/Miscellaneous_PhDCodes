%% Correction 

FileList = dir('L:\OT_ParticleWall\Unmodified\*z=1d*0.0100*.csv');
minsize = 1e5;
Vx = zeros(minsize,2*length(FileList));
for i=1:length(FileList)
    filename = [FileList(i).folder,'\',FileList(i).name];
    data = dlmread(filename);
    data(:,2) = round(data(:,2)*0.91,2);
    data(:,3) = round(data(:,3)*0.93,2);
    dlmwrite(['L:\OT_ParticleWall\',FileList(i).name],data,',');
end

