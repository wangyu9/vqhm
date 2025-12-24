function [t,l,h] = render_mesh(V,F,varargin)

% default values if not import from varargin
azel = [0 0];
C = [150,220,150]./255.*1.1;
w = [];
c_range = [];
uf = [];
nvar = length(varargin);
ii=1;
cmap = 'weights-neg';

LightSource = 'default';
FaceLighting = 'phong';
EdgeColor = 'none';
LightMultiplier = 1;
while(ii<=nvar)
   if(strcmp(varargin{ii},'view'))
       azel = varargin{ii+1};
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
   elseif(strcmp(varargin{ii},'ScaleColor'))
       w = varargin{ii+1};
       ii = ii + 1;
   elseif(strcmp(varargin{ii},'ColorMap'))
       cmap = varargin{ii+1};
       ii = ii + 1;    
   elseif(strcmp(varargin{ii},'ColorAxis'))
       c_range = varargin{ii+1};
       ii = ii + 1;
   elseif(strcmp(varargin{ii},'EdgeColor'))
       EdgeColor = varargin{ii+1};
       ii = ii + 1;    
   elseif(strcmp(varargin{ii},'LightMultiplier'))
       LightMultiplier = varargin{ii+1};
       ii = ii + 1;            
   end
   ii = ii + 1;
end


%%

if ~isempty(uf)
    uf = (uf-min(uf))/(max(uf)-min(uf));

    assert(min(uf)>=0);
    assert(max(uf)<=1);
    %C = value2color(u);
    ccc = colormap(my_colormap(cmap));

    C = squeeze(ind2rgb(floor(uf*size(ccc,1))+1,ccc));
    
    C = int8(255*C)+1; % this is a temp hack, do it right later. 
end


if(length(w)==0)
    t = tsurf(F,V,'EdgeColor',EdgeColor,'FaceColor',C','FaceLighting',FaceLighting,'SpecularStrength',1.);
else
    t = tsurf(F,V,'EdgeColor',EdgeColor,'FaceColor','interp','FaceLighting',FaceLighting);
    %t = trisurf(F,V(:,1),V(:,2),V(:,3));%,'EdgeColor',EdgeColor,'FaceColor','interp','FaceLighting',FaceLighting);
    %if strcmp(cmap,'weights-neg')
    colormap(my_colormap(cmap));
    %else
    %    colormap(cmap);
    %end
    %caxis([-0.2,1]);
    if(~isempty(c_range))
       caxis([c_range(1),c_range(2)]); 
    end
    %colorbar
    if(length(w)~=0)
        set(t,'CData',w);
    end
    drawnow;
end

l = [];
h = [];

% %set(t,'position',[100 100 640 360]);
% view(azel); % http://www.mathworks.com/help/matlab/ref/view.html
% 
% grid off;
% axis off;
% camzoom(1);
% 
% set(gcf,'color','w');
% %whitebg([1,1,1]);
% 
% %lightangle(-45,30)
% 
% 
% lightPos = [mean(V(:,1:2)), mean(V(:,3)) + mean((max(V(:,1:2))-min(V(:,1:2))))];
% 
% switch(LightSource)
%     case 'none'
%     case 'default'
% 
%         s = 0.32*LightMultiplier;
% 
%         l = [ 
%             light('Position',[+1 +1 +1],'Style','infinite','Color',s*[1,1,1]);
%             light('Position',[+1 +1 -1],'Style','infinite','Color',s*[1,1,1]);
%             light('Position',[+1 -1 +1],'Style','infinite','Color',s*[1,1,1]);
%             light('Position',[+1 -1 -1],'Style','infinite','Color',s*[1,1,1]);
%             light('Position',[-1 +1 +1],'Style','infinite','Color',s*[1,1,1]);
%             light('Position',[-1 +1 -1],'Style','infinite','Color',s*[1,1,1]);
%             light('Position',[-1 -1 +1],'Style','infinite','Color',s*[1,1,1]);
%             light('Position',[-1 -1 -1],'Style','infinite','Color',s*[1,1,1]);
%             %mean(V(:,1:2)), 4*(max(V(:,1:2))-min(V(:,1:2)))
%             %% light('Position',[10 -10  13],'Style','local')
%           ];
% end
%    
%     camproj('persp');
%     axis equal;
%     h = [];
%     %h = my_add_shadow(t,l,'Ground',[0,0,-1,-50]);
% %}