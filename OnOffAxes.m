function OnOffAxes(lim,ax,opt)
% OnOffAxes Adjust axes to display ON-OFF ISR cycles

arguments
  lim (:,1)
  ax = gca
  opt.ticks (1,:) = [0.3,0.6] % for polar
  opt.grid (1,1) {mustBeLogical} = true % for polar
  opt.repeat (1,1) {mustBeLogical} = true % for linear
end

% default value
if isscalar(lim)
  lim = [0,lim];
end

for i = 1 : numel(ax)

  if isa(ax(i),'matlab.graphics.axis.Axes')

    % prepare coordinates
    y_text = lim(1) - 0.075 * diff(lim);
    if opt.repeat
      x_max = 4*pi;
      x_text = [pi/2,5*pi/2;3*pi/2,7*pi/2];
      y_text = y_text * [1,1];
      x_line = [0,pi,NaN,2*pi,3*pi;pi,2*pi,NaN,3*pi,4*pi];
      y_line = lim(1) * ones(1,5);
    else
      x_max = 2*pi;
      x_text = [pi/2;3*pi/2];
      x_line = [0,pi;pi,2*pi];
      y_line = lim(1) * ones(1,2);
    end

    % adjust limits
    set(ax(i),'XLim',[0,x_max],'YLim',lim,'XColor','none')
    text(ax(i),x_text(1,:),y_text,'ON','Color',.4*ones(3,1),'FontSize',ax(i).FontSize,'FontWeight','bold','HorizontalAlignment','center','VerticalAlignment','top')
    text(ax(i),x_text(2,:),y_text,'OFF','Color','#c9a6f7','FontSize',ax(i).FontSize,'FontWeight','bold','HorizontalAlignment','center','VerticalAlignment','top')

    % plot axes
    h(1) = line(ax(i),x_line(1,:),y_line,'Color',.4*ones(3,1),'LineWidth',2);
    h(2) = line(ax(i),x_line(2,:),y_line,'Color','#c9a6f7','LineWidth',2);

  else

    % remove MATLB axes, adjust limits
    set(ax(i),'RLim',lim,'ThetaLim',[0,360],'RTick',[],'ThetaTick',[],'RGrid','off','ThetaGrid','off')

    % plot circles and cross
    if opt.grid
      h = matlab.graphics.chart.primitive.Line.empty;
      if ~isempty(opt.ticks)
        h(1:numel(opt.ticks)) = polarplot(ax(i),linspace(0,2*pi,100),opt.ticks.*ones(100,size(opt.ticks,2)),'Color',[0.8,0.8,0.8],'LineStyle','--');
      end
      h(end+1) = polarplot(ax(i),[0;pi;NaN;pi/2;3*pi/2],ones(5,1),'Color',[0.8,0.8,0.8]);
      for tick = opt.ticks
        text(ax(i),pi/2-0.1,tick*0.98,string(tick),'FontSize',ax(i).FontSize,'VerticalAlignment','bottom')
      end
    end

    % plot rotating arrows
    mksz = 100*ax(i).FontSize/15;
    h([numel(opt.ticks)+2,numel(opt.ticks)+3]) = polarplot(ax(i),linspace(-pi/2,pi/2,100),lim(2)*ones(100,2),'Color',.4*ones(3,1),'LineWidth',ax(i).LineWidth*2.8);
    h([end+1,end+2]) = polarplot(ax(i),linspace(pi/2,3*pi/2,100),lim(2)*ones(100,2),'Color','#c9a6f7','LineWidth',ax(i).LineWidth*2);
    h(end+1) = polarscatter(ax(i),pi/2,lim(2),mksz,'<','filled','MarkerFaceColor',.4*ones(3,1));
    h(end+1) = polarscatter(ax(i),3*pi/2,lim(2),mksz,'>','filled','MarkerFaceColor','#c9a6f7');

    % ON, OFF and R labels
    text(ax(i),pi/6,lim(2)*1.05,'ON','Color',.4*ones(3,1),'FontSize',ax(i).FontSize,'FontWeight','bold')
    text(ax(i),7*pi/6,lim(2)*1.05,'OFF','Color','#c9a6f7','FontSize',ax(i).FontSize,'FontWeight','bold','HorizontalAlignment','right')

  end

  % remove items from legend
  RemoveFromLegend(h)
  clearvars h

end