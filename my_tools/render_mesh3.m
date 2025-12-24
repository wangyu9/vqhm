function [t,l,h] = render_mesh3(V,F,varargin)

% default values if not import from varargin
azel = [];
C = []; % do not change this.
%C = [150,220,150]./255.*1.1;
u = [];
uf = [];
c_range = [];
nvar = length(varargin);

AO = [];
cmap = 'weights-neg';
Shading = 1.0;
VertexColor = [];
LightSource = 'default';
FaceLighting = 'gouraud';%'phong';
EdgeColor = 'none';
LightMultiplier = 1;

VP = [];
VN = [];

ii=1;
while(ii<=nvar)
   if(strcmp(varargin{ii},'view'))
       azel = varargin{ii+1};
       ii = ii + 1;
   elseif(strcmp(varargin{ii},'AmbientOcclusion'))
       AO = varargin{ii+1};
       ii = ii + 1;   
   elseif(strcmp(varargin{ii},'Shading'))
       Shading = varargin{ii+1};
       ii = ii + 1;
   elseif(strcmp(varargin{ii},'LightSource'))
       LightSource = varargin{ii+1};
       ii = ii + 1;
   elseif(strcmp(varargin{ii},'FaceColor'))
       C = varargin{ii+1};
       ii = ii + 1;
   elseif(strcmp(varargin{ii},'FaceScaleColor'))
       uf = varargin{ii+1};
       ii = ii + 1;
   elseif(strcmp(varargin{ii},'VertexColor'))
       VertexColor = varargin{ii+1};
       ii = ii + 1;
   elseif(strcmp(varargin{ii},'ScaleColor'))
       u = varargin{ii+1};
       ii = ii + 1;
   elseif(strcmp(varargin{ii},'ColorMap'))
       cmap = varargin{ii+1};
       ii = ii + 1;    
   elseif(strcmp(varargin{ii},'ColorAxis'))
       c_range = varargin{ii+1};
       ii = ii + 1;
   elseif(strcmp(varargin{ii},'Quiver'))
       Quiver = varargin{ii+1};
       assert(size(Quiver,2)==6);
       VP = Quiver(:,1:3);
       VN = Quiver(:,4:6);
       ii = ii + 1;   
   elseif(strcmp(varargin{ii},'EdgeColor'))
       EdgeColor = varargin{ii+1};
       ii = ii + 1;    
   elseif(strcmp(varargin{ii},'LightMultiplier'))
       LightMultiplier = varargin{ii+1};
       ii = ii + 1;
   else
       error(['Unknown parameter!\n']);
   end
   ii = ii + 1;
end

%%
if(size(V,2)==2)
   warning('z axis is missing! Append zero');
   V = [V,zeros([size(V,1),1])];
end
%%
if false
h = figure('Position',[100,450,800,900],'PaperPositionMode', 'auto');
%h = figure('PaperPositionMode', 'auto');
set(gcf,'papersize',[16 12])
end
% h = [];
%%

%t = tsurf(F,V,'EdgeColor',EdgeColor,'FaceColor','interp','FaceLighting',FaceLighting);
R = axisangle2matrix([0 0 1],pi)*axisangle2matrix([1 0 0],pi/2);
VD = V*R;

if isempty(C) && isempty(uf)
    t = trisurf(F,VD(:,1),VD(:,2),VD(:,3),'EdgeColor',EdgeColor,'FaceColor','interp','FaceLighting',FaceLighting);
    %t.Vertices = V*axisangle2matrix([0 0 1],pi)*axisangle2matrix([1 0 0],pi/2);

    if(~isempty(u))
        u = (u-min(u))/(max(u)-min(u));
        assert(min(u)>=0);
        assert(max(u)<=1);
        if 0 
            set(t,'CData',u);
            t.CDataMapping = 'direct';
            caxis([0,1])
        else      
            %VertexColor = value2color(u);
            VertexColor = squeeze(ind2rgb(floor(u*size(colormap,1))+1,colormap));
        end
    end

    if(~isempty(VertexColor))
        try
        if(isempty(AO))
            AO = ambient_occlusion(V,F,V,per_vertex_normals(V,F),1000);
        end
        catch
           AO = ones(size(V,1),1); 
        end

        %t.FaceVertexCData = bsxfun(@times,VertexColor,1-Shading*AO);
        tmp = bsxfun(@times,VertexColor,1-Shading*AO);
        set(t,'FaceColor','interp',...
           'CData',tmp,...
           'CDataMapping','scaled');
        %t.CData = u;%bsxfun(@times,VertexColor,1-Shading*AO);
    end

else
    t = trisurf(F,VD(:,1),VD(:,2),VD(:,3),'EdgeColor',EdgeColor,'FaceColor','flat','FaceLighting',FaceLighting);
    if isempty(uf)
        assert(size(C,1)==size(F,1));
        FaceColor = C;
    else
        uf = (uf-min(uf))/(max(uf)-min(uf));
        assert(min(uf)>=0);
        assert(max(uf)<=1);
        FaceColor = squeeze(ind2rgb(floor(uf*size(colormap,1))+1,colormap));
    end
    %t = trisurf(F,VD(:,1),VD(:,2),VD(:,3),'EdgeColor',EdgeColor,'FaceColor','flat','FaceLighting',FaceLighting,'CData',C');
    set(t,'FaceColor','flat',...
       'FaceVertexCData',FaceColor,...
       'CDataMapping','scaled');
end

colormap(my_colormap(cmap));
% if(~isempty(c_range))
%    caxis([c_range(1),c_range(2)]); 
% end

t.SpecularStrength = .4;%0.3;
t.DiffuseStrength = .45;%0.1;
t.AmbientStrength = .6;%0.7;
t.SpecularColorReflectance = .3;
t.SpecularExponent = 7;

switch(LightSource)
    case 'none'
    case 'default'

        s = 0.3*LightMultiplier;

        l = [ 
            light('Position',[+1 +1 +1],'Style','infinite','Color',s*[1,1,1]);
            light('Position',[+1 +1 -1],'Style','infinite','Color',s*[1,1,1]);
            light('Position',[+1 -1 +1],'Style','infinite','Color',s*[1,1,1]);
            light('Position',[+1 -1 -1],'Style','infinite','Color',s*[1,1,1]);
            light('Position',[-1 +1 +1],'Style','infinite','Color',s*[1,1,1]);
            light('Position',[-1 +1 -1],'Style','infinite','Color',s*[1,1,1]);
            light('Position',[-1 -1 +1],'Style','infinite','Color',s*[1,1,1]);
            light('Position',[-1 -1 -1],'Style','infinite','Color',s*[1,1,1]);
          ];
end
   

    %h = my_add_shadow(t,l,'Ground',[0,0,-1,-50]);




camproj('persp');
axis equal;
grid off;
if(~isempty(azel))
view(azel); % http://www.mathworks.com/help/matlab/ref/view.html
end
axis off;
%camzoom(1);

set(gcf,'color','w');


if(~isempty(VP))
   VP = VP * R;
   VN = VN * R;
   quiver3(VP(:,1),VP(:,2),VP(:,3),VN(:,1),VN(:,2),VN(:,3),'ShowArrowHead','off'); 
   % does not work. quiver3(VP(:,1)-0.5*VN(:,1),VP(:,2)-0.5*VN(:,2),VP(:,3)-0.5*VN(:,3),VN(:,1),VN(:,2),VN(:,3),'ShowArrowHead','off'); 
end


if (max(V(:,3))-min(V(:,3))<1e-10)
    view([0,0]);
end

drawnow;
