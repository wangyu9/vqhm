function [table_results] = json_to_table(JS)

fn = fieldnames(JS);

%%
variable_type = ["string",   "double",          "int64"]; 
variable_name = ["name",     "time",   "flips"]; 

%%

table_results = table('Size',[numel(fn),numel(variable_type)],'VariableTypes',variable_type,'VariableNames',variable_name);

%%

for k=1:numel(fn)

    name = fn{k};

    table_results(k, 'time') = {sum(JS.(name).time)}; 

    table_results(k, 'name') = {name}; 

    table_results(k, 'flips') = {JS.(name).nInverted(end)}; 

end

