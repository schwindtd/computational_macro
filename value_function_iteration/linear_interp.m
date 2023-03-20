function [w] = linear_interp(x,y,z)
%Linear interpolate on grid
%    Given original x vector, y vector of new points and z an original 
%    matrix, computes a new matrix w, with size (length(y), ncol(z))
s = size(z);
w = zeros(length(y), s(2));
for i = 1:s(2)
    for j=1:length(y)
        if (j==1)
           w(j,i) = z(j,i);
        elseif (j==length(y))
            w(j,i) = z(s(1),i);
        else
            leftn = x(y(j) > x);
            rightn = x(y(j) < x);
            ml = max(leftn);
            mr = min(rightn);
            left = find(x == ml);
            right = find(x == mr);
            w(j, i) = z(left, i) + ...
            (z(right, i) - z(left, i))/...
            (x(right) - x(left))...
            *(y(j) - x(left));
        end
    end
end
           
end

