function [] = draw_uv_with_distortion(u,v,F,uf)

f = size(F,1);

% flipped = doublearea([u,v], F) < 0;
% 
% 
% C = ([1,1,1].*ones(f,3)) .* (1-flipped(:)) + ([0.8,0.2,0.2].*ones(f,3)) .* flipped(:);
% [t,l,h] = render_mesh2([u,v],F,'EdgeColor',[0,0,0],'FaceColor',[1,1,1]);
% set(t,'LineWidth',0.9);
% set(t,'FaceVertexCData',C);
% axis off;
% axis equal;
% hold on

%%
[t,l,h] = render_mesh2([u,v],F,'EdgeColor',[0,0,0]); %,'ColorMap','heat','FaceScaleColor',uf);
%%
C = uf.*[1,0,0]+[1,1,1].*(1-uf);


set(t,'FaceColor','flat',...
       'FaceVertexCData',C,...
       'CDataMapping','scaled');

% set(t,'LineWidth',0.9);
% set(t,'FaceVertexCData',C);
axis off;
axis equal;

% [t,l,h] = render_mesh2([u,v],F(flipped,:),'EdgeColor',[0,0,0],'FaceColor',[0.8,0.2,0.2], 'ColorMap','heat','FaceScaleColor',uf);


% [t,l] = render_mesh3([u,v],F,'EdgeColor',[0,0,0],'ColorMap','heat','FaceColor',uf.*[1,0,0]+[1,1,1].*(1-uf));

% [t,l] = render_mesh3([u,v],F,'EdgeColor',[0,0,0],'ColorMap','heat','FaceScaleColor',uf);