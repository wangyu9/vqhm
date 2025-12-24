
%%

for ii=1:numel(history)
%%
    u = history{ii}.u;
    v = history{ii}.v;
    
    [t,l,h] = render_mesh2([u,v],F,'EdgeColor',[0,0,0],'FaceColor',[1,1,1]); 
    axis off; 
    axis equal; 
    
    pause(2);
    
    legend(['frame: ',num2str(ii)]);


end