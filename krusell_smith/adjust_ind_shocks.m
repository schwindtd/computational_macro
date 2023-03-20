function [indshocks, indshocks_idx] = adjust_ind_shocks(indshocks,aggshocks_idx,u_g, u_b)
dim = size(indshocks);
% Adjust shocks so aggregate unemployment rates exactly match u_g and u_b
unemp_real = 1- mean(indshocks,2);
unemp_true = (aggshocks_idx==1)*u_b + (aggshocks_idx==2)*u_g;
unemp_diff = unemp_true - unemp_real;
num2change = floor(abs(unemp_diff)*dim(2));
dir2change = sign(unemp_diff);

for i=1:dim(1)
    if num2change(i) > 0
        if dir2change(i) == 1 % Need to flip some individuals to unemployed
            emp_idx = find(indshocks(i,:) == 1);
            rand_idx = randperm(length(emp_idx),num2change(i));
            flip_idx = emp_idx(rand_idx);
            indshocks(i, flip_idx) = 0;
        elseif dir2change(i) == -1 % Need to flip some individuals to employed
            emp_idx = find(indshocks(i,:) == 0);
            rand_idx = randperm(length(emp_idx),num2change(i));
            flip_idx = emp_idx(rand_idx);
            indshocks(i, flip_idx) = 1;
        end
    end
end

indshocks_idx = indshocks + 1;
end

