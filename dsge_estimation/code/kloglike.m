function ll = kloglike(F,v)
%Function computing log likelihood using outputs from Kalman Filter
dim_v = size(v);
n = dim_v(1);
p = dim_v(2);
c = -(n*p/2)*log(2*pi);
s=0;
for i=1:n
s = s + log(det(F(:,:,i))) + v(i,:)*(F(:,:,i)^-1)*v(i,:)';
end
ll = c - 0.5*s;

end

