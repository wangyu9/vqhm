function [p] = bold_hist2(stats,varargin)
xscale = 'linear';
%%
p = [];

h = figure('position',[0,0,800,600]);


set(gcf,'papersize',[16 12])


%%

hist = histogram(stats, 'BinWidth', 0.001); 

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
