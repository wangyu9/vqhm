function [names] = find_files_with_ext(folder)

files = dir(folder);
names = {files(3:end).name}; 
isMesh = false(numel(names),1);

for k=1:numel(names)
    [filepath,name,ext] = fileparts(names{k});
    if(strcmp(ext,'.obj'))
        isMesh(k) = true;
    end
end
names = names(isMesh);