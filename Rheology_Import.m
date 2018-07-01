close all; clear variables

[FileName,PathName] = uigetfile('*.csv');
fileID = fopen([PathName,FileName],'r');
Intro = textscan(fileID,'%s',18,'Delimiter','\n'); %Read 18 header lines


Block  = 1;

while (~feof(fileID))                               % For each block:                         
   InputText = textscan(fileID,'%s',12,'delimiter','\n');  % Read 12 header lines
   HeaderLines{Block,1} = InputText{1};
   fprintf('Block: %s\n', num2str(Block))           % Print block number to the screen
   FormatString = [repmat('%q',1,6),'%*[^\n]'];     % Create format string
   InputText = textscan(fileID,FormatString,20, ...    % Read data block
      'delimiter',',');
   
   InputText = [InputText{2},InputText{5},InputText{6}]; %Select relevant columns and concatenate them in an array
   Data{Block,1} = cellfun(@str2num, InputText);              
   [NumRows,NumCols] = size(Data{Block});           % Determine size of table
   disp(cellstr(['Table data size: ' ...
      num2str(NumRows) ' x ' num2str(NumCols)]));
   disp(' ');                                       % New line
   
   eob = textscan(fileID,'%s',1,'delimiter','\n');  % Read and discard end-of-block marker 
   Block = Block+1;                                 % Increment block index
end

fclose(fileID);
%% 
figure(1)
loglog(Data{1,1}(:,1)/100,Data{1,1}(:,2),'or','MarkerFaceColor','r')
hold on
loglog(Data{1,1}(:,1)/100,Data{1,1}(:,3),'or')
loglog(Data{2,1}(:,1)/100,Data{2,1}(:,2),'ob','MarkerFaceColor','b')
loglog(Data{2,1}(:,1)/100,Data{2,1}(:,3),'ob')
hold off

% figure(2)
% loglog(100*Data{3,1}(:,1),Data{3,1}(:,2),'or','MarkerFaceColor','r')
% hold on
% loglog(100*Data{3,1}(:,1),Data{3,1}(:,3),'or')
% loglog(100*Data{4,1}(:,1),Data{4,1}(:,2),'ob','MarkerFaceColor','b')
% loglog(100*Data{4,1}(:,1),Data{4,1}(:,3),'ob')
% hold off

