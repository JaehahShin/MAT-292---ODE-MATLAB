function [t, y] = lab3_shinja19_adv_eul(f, t0, tN, y0, h)
    t = [t0];
    y = [y0];
    tol = 1e-8;
    function [Y, Z, D] = helper(y, h, t) % Helperfunction with imp_eulder 
        Y = y + h * f(t, y);
        S_L = f(t, y);
        S_R = f(t + h, Y);
        improv_slope = (S_L + S_R) / 2;
        Y = y + h * improv_slope;
 
        Z1_first = y + 0.5 * h * f(t, y);
        S_L = f(t, y);
        S_R = f(t + 0.5 * h, Z1_first);
        improved_slope = (S_L + S_R) / 2;
        Z = Z1_first + 0.5 * h * improved_slope;
 
        S_L = f(t + 0.5 * h, Z1_first);
        S_R = f(t + h, Z1_first);
        improved_slope = (S_L + S_R) / 2;
        Z = Z1_first + 0.5 * h * improved_slope;
        D = Z - Y;
    end
    
    while t(end) < tN
        [Y, Z, D] = helper(y(end), h, t(end));
        while abs(D) >= tol % check error range 
            h = 0.9 * h * min(max(tol / abs(D), 0.3), 2);
            [Y, Z, D] = helper(y(end), h, t(end));
        end
        % Now everything is good 
        y = [y, Z + D];
        t = [t, t(end) + h];
    end
end