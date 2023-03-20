function r = right(k, vect2)
    if (k < min(vect2))
        r=2;
    elseif (k > max(vect2))
        r=length(vect2);
    else
        r = find(vect2 == min(vect2(vect2>k)));
    end
end
