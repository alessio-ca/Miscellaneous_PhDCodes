function [] = csv_data_structuring(filename)

% This program takes an original track.csv files as produced by Analyse (on
% the POM 0.26 microscope) and restructures it to be used for the power
% spectrum microrheology routines.
% It assumes a N-by-3 input file, with format:
% 1   - time
% 2:3 - X and Y positions (in pixels)
%
% INPUTS
%
% filename - The filename of the track file.
%
% OUTPUTS
%
% Creates restructured track file and save it in "Trackfile.mat", suitable for
% input in the microrheo routines.
%
% Output format:
%
% 1:2 - X and Y positions (in micrometers)
% 3   - frame # (conversion to time is executed later)
% 4   - Bead ID (trivially 1, but needed for consistency)


res = csvread(filename);
res = circshift(res,[0 -1]);
res(:,3) = 1:size(res,1);
res(:,4) = ones(size(res,1),1);
save('Trackfile.mat', 'res' );



end
