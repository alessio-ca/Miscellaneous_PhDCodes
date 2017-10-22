function C = DWSread(filename,type)
% type = 'raw' or 'MR'
% created by Zhongyang Xing 21/05/2017

fileID = fopen(filename);

switch type
    case 'raw'
        C = textscan(fileID,'%f %f','HeaderLines',21);
        % column1: Time lag 
        % column2: g2-1
    case 'MR'
        C = textscan(fileID,'%f %f %f %f %f %f %f %f %f %f','HeaderLines',21);
        % column1: Lag Time (s)
        % column2: MSD (micron^2)
        % column3: Fit Lag Time (s)
        % column4: MSD Fit (micron^2)
        % column5: Frequency (rad/s)
        % column6: G1(Pa)
        % column7: G2(Pa)
        % column8: G*(Pa)
        % column9: Complex Viscosity (Pa s)
        % Loss Tangent (-)
end
C = cell2mat(C);
fclose(fileID);
end 