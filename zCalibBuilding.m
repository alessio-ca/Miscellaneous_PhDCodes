% direcZseries='F:\AC_2016_05_30\';
% mkdir(direcZseries,'dircalib');
% listSeries=dir([direcZseries,'\zcalib\','60*']);
% 
% for i=1:length(listSeries)
%     copyfile([direcZseries,'\zcalib\',listSeries(i).name,'\*_000000.tiff'],[direcZseries,'\dircalib'])
% end

direcZseries='F:\AC_2017_05_24\';
picdir = 'ME_SiOil_F108_12VStir\';
framedir = 'ME_F108_SiOil_Frames';
mkdir(direcZseries,framedir);
listSeries=dir([direcZseries,picdir,'10*']);
listSeries=listSeries([listSeries.isdir]);

for i=1:length(listSeries)
    copyfile([direcZseries,picdir,listSeries(i).name,'\*_000000.tiff'],[direcZseries,framedir])
end
