function W = hat(w)
    % This function implements the ^-operator.
    % It converts a 3x1-vector into a skew symmetric matrix.
    [r, c] = size(w);
    if r == 3 && c == 1
        W = zeros(3);
        W(1,2) = -w(3);
        W(2,1) = w(3);
        W(1,3) = w(2);
        W(3,1) = -w(2);
        W(2,3) = -w(1);
        W(3,2) = w(1);
    else
        error('Variable w has to be a 3-component vector!')
    end
end