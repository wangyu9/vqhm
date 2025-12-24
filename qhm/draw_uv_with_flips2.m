function [] = draw_uv_with_flips2(u,v,F,args)

f = size(F,1);

flipped = doublearea([u,v], F) < 0;

% flipped = V(F(:,1),1) > 0;
C = ([1,1,1].*ones(f,3)) .* (1-flipped(:)) + ([0.8,0.2,0.2].*ones(f,3)) .* flipped(:);
% C = C * 255;
[t,l,h] = render_mesh2([u,v],F,'EdgeColor',[0,0,0],'FaceColor',[1,1,1]);
set(t,'LineWidth',0.9);
set(t,'FaceVertexCData',C);
axis off;
axis equal;
hold on

[t,l,h] = render_mesh2([u,v],F(flipped,:),'EdgeColor',[0,0,0],'FaceColor',[0.8,0.2,0.2]);