function [deriv] = numderiv(grid, vect)
%Computes numerical derivatives over a grid
%   Derivative = (Left Derivative + Right Derivative)/2, except if at the
%   edges where instead the derivative is simply the slope from the edge to
%   the nearest non-edge point
s = size(grid);
for i = 1:s(2)
    for j = 1:s(1)
        if (j == 1)
            num = grid(j+1, i) - grid(j, i);
            denom = vect(j+1) - vect(j);
            right_deriv = num/denom;
            deriv(j, i) = right_deriv;
        elseif (j == s(1))
            num = grid(j, i) - grid(j-1, i);
            denom = vect(j) - vect(j-1);
            left_deriv = num/denom;
            deriv(j, i) = left_deriv;
        else
            lnum = grid(j, i) - grid(j-1, i);
            ldenom = vect(j) - vect(j-1);
            rnum = grid(j+1, i) - grid(j, i);
            rdenom = vect(j+1) - vect(j);
            left_deriv = lnum/ldenom;
            right_deriv = rnum/rdenom;
            deriv(j, i) = (left_deriv + right_deriv)/2;
        end
        
    end
end
end

