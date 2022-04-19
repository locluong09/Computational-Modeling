
function out = linearInterpolation(x, y, z)
    %{
    This function is not the same with matlat interp1, it just return the
    value of new value of fluid properties when get the value of z
    z is a new value of pressure that need to be interpolated to find the values
    of fluid properties
    %}
    n = length(x);
    if ismember(z,x)
        out = y(find(x == z));
        return
    end
    % finding the interval that contains the value of z
    left = 1;
    right = n;
    while left < right
        if z < (x(left) + x(right)) / 2
            right = right -1;
        else
            left = left + 1;
        end
    end
    % update left and right because after the while loop right = left
    if z > x(right)
        right = right +1;
    else
        left = left - 1;
    end
    
    % using linear interpolation to update the value of z
    out = (z - x(left))*(y(right) - y(left))/(x(right) - x(left)) + y(left);
end


