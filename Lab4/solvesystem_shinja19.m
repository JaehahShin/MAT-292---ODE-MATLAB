function [time, x] = solvesystem_shinja19(f, g, t0,tN,x0,h)
%  SOLVESYSTEM_SHINJA19 x1'=f(t,x1,x2), x2'=g(t,x1,x2) where t0 and tN are the start and end points 
% of the interval on which to solve the ODE, h is the stepsize, and x0 is a vector 
% for the initial condition of the system of ODEs x(t0)=x0 return a row vector 
% of times and a matrix of approximate solution values (the first row has the 
% approximation for x1 and the second row has the approximation for x2). Note: 
% you will need to use a loop to do this exercise. You will also need to recall 
% the Heun/Improved Euler algorithm learned in lectures.
 len = round((tN-t0)/h,0);
 time = linspace(t0, tN, len);
 x = zeros(2, len);
 % initial Vlues:
 time(1) = t0;
 x(:,1) = x0; % at(:,1) % gives the first column (and so on) which will be x0
    for i = 2:len
        time(i) = time(i-1) + (i-1)*h;
        % First
        x(1,i) = x(1, i-1) + h*f(time(i-1), x(1,i-1), x(2,i-1));
        x(2,i) = x(2,i-1) + h*g(time(i-1), x(1,i-1), x(2,i-1));
        % slope of left
        S_L = [f(time(i-1), x(1,i-1), x(2,i-1));
                g(time(i-1), x(1,i-1),x(2,i-1))];
        % right one
        S_R = [f(time(i), x(1,i), x(2,i));
                g(time(i),x(1,i),x(2,i))];
        improv_slope = (S_L + S_R) ./ 2;
        x(:,i) = x(:, i -1) + h*improv_slope;
        
    end 
end