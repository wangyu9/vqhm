function [val] = json_reader(fname)

%%
if false
%%
JS = json_reader('50lines_records.json');
%%

JS.lucy_o_G.time ./ JS.lucy_o_G.LBFGSIternum
end



%  fname = 'meta_Electronics.json'; 

fid = fopen(fname); 
raw = fread(fid,inf); 
str = char(raw'); 
fclose(fid); 
val = jsondecode(str);