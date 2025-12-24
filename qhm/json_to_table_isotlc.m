function [table_results] = json_to_table_isotlc(JS)

fn = fieldnames(JS);

%%
variable_type = ["string",   "double",          "int64"]; 
variable_name = ["name",     "garanzha_time",   "garanzha_flips"]; 

%%

table_results = table('Size',[numel(fn),numel(variable_type)],'VariableTypes',variable_type,'VariableNames',variable_name);

%%

for k=1:numel(fn)

    name = fn{k};

    table_results(k, 'garanzha_time') = {sum(JS.(name).time)}; 

    table_results(k, 'name') = {name}; 

    table_results(k, 'garanzha_flips') = {JS.(name).nInverted(end)}; 

end

