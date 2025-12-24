function [pathname,fullname,id] = find_max_with_prefix(prefix,folder,ext)

len = numel(prefix); 

len_ext = numel(ext); 


% startIndex = regexp('result_Du_orgiter_95','result_Du_orgiter_\d*')

files = dir(folder); 

id = -1;

for i=1:numel(files)
    name = files(i).name;
    if numel(name) > (len+len_ext) && strcmp(name(1:len), prefix)
        current = str2num( name(len+1:end-len_ext) ); 
        id = max(id, current); 
    end
end



fullname = [prefix, num2str(id), ext]; 

pathname = folder + '/' + fullname;