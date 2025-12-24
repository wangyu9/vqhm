function [p] = plot_stats(stats,varargin)

metric = [];
xscale = 'log';
params_to_variables = containers.Map( { 'metric', 'xscale'},...
                                       {'metric', 'xscale'});


%% Shared Parsing Code Segment

v = 1;
while v <= numel(varargin)
    param_name = varargin{v};
    if isKey(params_to_variables,param_name)
        assert(v+1<=numel(varargin));
        v = v+1;
        % Trick: use feval on anonymous function to use assignin to this workspace 
        feval(@()assignin('caller',params_to_variables(param_name),varargin{v}));
    else
      error('Unsupported parameter: %s',varargin{v});
    end
    v=v+1;
end
%%
p = [];

h = figure('position',[0,0,800,600]);


set(gcf,'papersize',[16 12])


%%

hist = histogram(stats.sigma_ratio); 

if ~isempty(xscale)
set(gca, 'xscale',xscale);
end

%%

set(gca,'fontsize',20,'linewidth',2)

box on;

% make the border tight to save space
% https://www.mathworks.com/matlabcentral/newsreader/view_thread/156326
if 1
    % tight
    set(gca,'LooseInset',get(gca,'TightInset'))
end
% https://www.mathworks.com/matlabcentral/answers/12987-how-to-save-a-matlab-graphic-in-a-right-size-pdf
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

p.h = h;
p.hist = hist;

%%

if false
p = []

h = figure('Position',[0,0,1440,1080],'PaperPositionMode', 'auto');
%h = figure('PaperPositionMode', 'auto');
set(gcf,'papersize',[16 12])

set(h,'defaultAxesColorOrder',[get(h,'defaultAxesColorOrder');0.25,0.25,0.25])

if ~isempty(line_colors)
   set(h,'defaultAxesColorOrder',line_colors);
end

set(gca,'fontsize',20,'linewidth',2)

hold on;

if size(u,1)==1
   u = u'; 
end

c = size(u,2);
ls = cell(c,1);
for i=1:c
    ls{i} = plot(0:size(u,1)-1,u(:,i),'LineWidth',3);
end

%ylim([0,35])

%h_legend = legend('','');
%set(h_legend,'FontSize',25,'Location','best');

box on;

% make the border tight to save space
% https://www.mathworks.com/matlabcentral/newsreader/view_thread/156326
if 0
    % tight
    set(gca,'LooseInset',get(gca,'TightInset'))
end
% https://www.mathworks.com/matlabcentral/answers/12987-how-to-save-a-matlab-graphic-in-a-right-size-pdf
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

p.ls = ls;
p.h = h;

end