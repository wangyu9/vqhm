
function [] = print_frames(history, F, frames, path)


for ii=1:numel(frames)
%%
    index = frames(ii);
    
    u = history{index}.u;
    v = history{index}.v;
    
    
    figure('position',[0 0 1000 1000]);
    [t,l,h] = render_mesh2([u,v],F,'EdgeColor',[0,0,0],'FaceColor',[1,1,1]);
    set(t,'LineWidth',2);
    axis off;
    axis equal;

    fname = [path,'/',num2str(index)];

    savefig([fname,'.fig'])
    print(fname,'-dpng');
    %im = imread(fname);
    %im = imrotate(im,-90);
    %imwrite(im,fname);
    %image_white2none([fname,'.png'],[fname,'.png']);
    close();
        

end