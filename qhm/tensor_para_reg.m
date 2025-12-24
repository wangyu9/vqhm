function [reg_value,reg_grad, rau_value, rau_grad] = tensor_para_reg(FA,dim,para_type)

% regularizer in at
reg_value = @(at) 0;
reg_grad  = @(at) 0;

% regularizer in au
rau_value = @(au) 0;
rau_grad  = @(au) 0;

Area = FA;