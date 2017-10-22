function [Dy,y] = DI(x)
%DI Dispersity Index calculator for a set of data

Dy = std(x)/mean(x);
y = mean(x);


