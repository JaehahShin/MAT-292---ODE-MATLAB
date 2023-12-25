% Use for loop to do this. 
% for loop

function [time,Y] = lab3_shinja19_imp_euler(f, t0, tN, y0, h)
    len = round((tN-t0)/h,0);
    time = linspace(t0, tN, len);
    Y = zeros(1, len);
    Y(1) = y0;

    for i = 2:len
        time(i) = t0 + (i-1)*h;
        Y(i) = Y(i-1) + h*f(time(i-1), Y(i-1));
        S_L = f(time(i-1), Y(i-1));
        S_R = f(time(i), Y(i));
        improv_slope = (S_L + S_R) / 2;
        Y(i) = Y(i -1) + h*improv_slope;
    end 

    %plot(time, Y)
    %title("Exercies 1: ImprovedODE Solver");
    %xlabel("t");
    %ylabel("y");
    
end