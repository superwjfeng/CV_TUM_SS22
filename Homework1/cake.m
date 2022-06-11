function Cake = cake(min_dist)
    % The cake function creates a "cake matrix" that contains a circular set-up of zeros
    % and fills the rest of the matrix with ones. 
    % This function can be used to eliminate all potential features around a stronger feature
    % that don't meet the minimal distance to this respective feature.
    
    range = -min_dist:min_dist;
    [x, y] = meshgrid(range, range);

    circle = sqrt(x.^2 + y.^2); 

    % returns a matrix containing a circular set-up of zeros and fills up
    % the rest of the matrix with ones
    Cake = ones(length(range));
    Cake(circle<=min_dist) = 0;
    Cake = logical(Cake);
end
