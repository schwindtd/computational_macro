function [shocks_z,shocks_z_idx] = agg_shocks_gen(T,z,P_z)
%% Aggregate Shocks Simulation
% Author: Daniel Schwindt
% Date: 12/28/2022

%% Generate Aggregate Productivity Chain
shocks_z = nan(T,1);
shocks_z(1) = z(1);
for i=2:T
    c_rand = rand;
    if (shocks_z(i-1) == z(1) && c_rand > P_z(1,2))
        shocks_z(i) = z(1);
    elseif (shocks_z(i-1) == z(1) && c_rand <= P_z(1,2))
        shocks_z(i) = z(2);
    elseif (shocks_z(i-1) == z(2) && c_rand > P_z(2,1))
        shocks_z(i) = z(2);
    else
        shocks_z(i) = z(1);
    end
end
% Shock indices
shocks_z_idx = (shocks_z==z(1)).*1 + (shocks_z==z(2)).*2;
end

