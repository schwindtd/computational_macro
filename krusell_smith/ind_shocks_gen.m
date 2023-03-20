function [shocks_e,shocks_e_idx] = ind_shocks_gen(T,I,u,P_z, P_ze, shocks_z_idx)
%% Generate Idiosyncratic Shock chain
shocks_e = nan(T,I);
c_rand = rand(T,I);
shocks_e(1,:) = (c_rand(1,:)>=u);
for i=2:T
    for j=1:I
        shocks_e(i,j) = (c_rand(i,j)<= ...
            P_ze(shocks_z_idx(i-1), ((shocks_e(i-1,j))+1),shocks_z_idx(i), 2));
    end
end
shocks_e_idx = shocks_e+1;
end

