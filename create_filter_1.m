function [filterFin]=create_filter_1(rMax,rMin)

filter1=1+mat2gray(fspecial('disk',rMax));
filter2=2*mat2gray(padarray(fspecial('disk',rMin),[rMax-rMin rMax-rMin],'both'));
filterFin=mat2gray(filter1-filter2);