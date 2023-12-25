function [time, y] = DE2_shinja19(f,t0,tN,y0,y1,h)
len = round(((tN-t0)/h));

time = t0:h:tN;
num_Points = length(time);
time_linspace = linspace(t0, tN, num_Points);

y = zeros(1,len);
y(1) = y0;
y(2) = y0+y1*h;

for i = 1:len-1
    % THis is first derivative of y
    y_one_p = (y(i+1) - y(i))/ h;
    % This is second derivative of y
    % calling function f with the current time, the first deriviateve which
    % is defined previously, and the current value of y. 
    y_doudble_p = f(time_linspace(i+2), y_one_p, y(i+1));
    y(i+2) = (h^2)*y_doudble_p + 2*y(i+1) - y(i);
end