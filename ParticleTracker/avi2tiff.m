function avi2tiff
[selectedFiles,PathName] = uigetfile('*.avi','Select the avi file(s): ','MultiSelect','on');
if isempty(selectedFiles)
    error('No input files.')
end

%If single entry, convert to cell
if ~iscell(selectedFiles)
    selectedFiles=cellstr(selectedFiles);
end
for i=1:length(selectedFiles)
   avi2tiff_func([PathName,selectedFiles{i}]);
end

function avi2tiff_func(file)

%Open avi file
video=VideoReader(file);
file=dir(file);
filename=file.name;
filefolder=file.folder;

%Create directory
myfolder=regexprep(filename,'.avi','');   
mkdir([filefolder,'/',myfolder]);
disp(myfolder);

%Print video info
fr=video.FrameRate;
disp(['Frame rate: ',num2str(fr)]);
fn=video.Duration*fr;
disp(['Number of frames: ',num2str(fn)]);

%Get frames
ii = 0;
while hasFrame(video)
   img = readFrame(video);
   imagename=[filefolder,'/',myfolder,'/',myfolder,'_',num2str(ii,'%.6d'),'.tiff'];
   imwrite(img,imagename)    
   ii = ii+1;
end
disp(' ')





