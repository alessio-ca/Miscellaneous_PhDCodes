h = gcf; %current figure handle
all_axesObjs = findobj(h, 'type', 'axes');
dataObjs = get(all_axesObjs, 'Children'); %handles to low-level graphics objects in axes
objTypes = get(dataObjs, 'Type');  %type of low-level graphics object
lineObjs = findobj(dataObjs, 'type', 'line');
xdata = get(dataObjs, 'XData');  %data from low-level grahics objects
ydata = get(dataObjs, 'YData');  %data from low-level grahics objects