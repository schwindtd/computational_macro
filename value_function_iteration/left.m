function l = left(k, vect2)
    if (k < min(vect2))
        l=1;
    elseif (k > max(vect2))
        l=length(vect2)-1;
    else
        l = find(vect2 == max(vect2(vect2<k)));
    end
end
