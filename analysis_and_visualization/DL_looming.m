% clearvars -except plotTable dataIDlist
datarefs = [2 3 1];
saveFig = 1;
hfig = figure;
backC = [1 1 1];
set(hfig,'position',[2062 264 866 633])
figPos = get(hfig,'position');
set(hfig,'color',backC,'position',figPos)
haxA = axes;
haxPos = get(haxA,'position');
growShiftAx = [-0.5,-0.3,0,0.15];%grow width, grow height, shift X, shift Y
haxPos = [haxPos(1)-haxPos(3)*growShiftAx(1)/2,haxPos(2)-haxPos(4)*growShiftAx(2)/2,...
    haxPos(3)*(1+growShiftAx(1)),haxPos(4)*(1+growShiftAx(2))]+[growShiftAx(3:4) 0 0];
xtickpos = [0 1 2 3];
yMin = -2.5;    yMax = 2.5;    yMid = 2.5;
ytickpos = (yMin:yMid:yMax);
ydrop = -0.1;
fontC = [0 0 0];
set(haxA,'xticklabel',[],'position',haxPos,'color','none','box','off','Xlim',[xtickpos(1)-0.5 xtickpos(end)+0.5],...
    'ticklength',[0.02 0.01],'YLim',[ytickpos(1) ytickpos(end)],'tickdir','out','xtick',xtickpos,...
    'yticklabel',[],'ytick',ytickpos,'ycolor','w','xcolor','w','nextplot','add')
plot([0 0],[ytickpos(1) ytickpos(end)],'color',[0 0 0],'linewidth',2,'parent',haxA)
xplot = cat(1,xtickpos,xtickpos,NaN(1,numel(xtickpos)));
yplot = cat(1,zeros(1,numel(xtickpos))+ydrop/2,zeros(1,numel(xtickpos))+ydrop,zeros(1,numel(xtickpos)));
plot([xtickpos(1) xtickpos(end)],[ydrop ydrop]/2,'k','linewidth',2,'parent',haxA(1))
colors = {'k','k','r'};
for iterL = 1:numel(datarefs)
    ydata = plotTable.returnData{dataIDlist{datarefs(iterL)}};
    hbox = boxplot(haxA,ydata,'positions',iterL-0.5,'boxstyle','outline',...
        'colors',colors{datarefs(iterL)},'notch','on','widths',0.75,'orientation','vertical');
    set(hbox(7,:),'marker','.','markeredgecolor','k','markersize',12)
    obj = findobj(hfig,'Tag','Upper Adjacent Value');
    set(obj,'visible','off')
    obj = findobj(hfig,'Tag','Lower Adjacent Value');
    set(obj,'visible','off')
    obj = findobj(hfig,'Tag','Upper Whisker');
    set(obj,'LineStyle','-','LineWidth',2)%,'color','k');
    obj = findobj(hfig,'Tag','Lower Whisker');
    set(obj,'LineStyle','-','LineWidth',2)%,'color','k');
    obj = findobj(hfig,'Tag','Box');
    set(obj,'LineWidth',2)
    obj = findobj(hfig,'Tag','Median');
    set(obj,'LineWidth',2)
    text(xtickpos(iterL+1)-0.75,yMin-2,{'median: ',num2str(median(ydata),2)},'rotation',0,...
        'color',fontC,'horizontalalignment','left','fontsize',14,'interpreter','none');
    text(xtickpos(iterL+1)-0.75,yMin-3,['n = ' num2str(numel(ydata))],'rotation',0,...
        'color',fontC,'horizontalalignment','left','fontsize',14,'interpreter','none');
end
set(haxA,'xticklabel',[],'position',haxPos,'color','none','box','off','Xlim',[xtickpos(1)-0.5 xtickpos(end)+0.5],...
        'ticklength',[0.02 0.01],'YLim',[ytickpos(1) ytickpos(end)],'tickdir','out','xtick',xtickpos,...
        'yticklabel',[],'ytick',ytickpos,'ycolor','w','xcolor','w','nextplot','add')
    
for iterL = 1:numel(datarefs)
    text(xtickpos(iterL+1)-0.75,yMin-1,(plotTable.plotID{dataIDlist{datarefs(iterL)}}),'rotation',0,...
        'color',fontC,'horizontalalignment','left','fontsize',14,'interpreter','none');
end
for iterL = 1:numel(ytickpos)
    text(xtickpos(1)-(0.1)*range(xtickpos),ytickpos(iterL),num2str(ytickpos(iterL)),'rotation',0,...
        'color',fontC,'horizontalalignment','right','fontsize',20,'interpreter','none');
end

ytickPOS = (ytickpos);
for iterT = 1:numel(ytickPOS)
    text(-0.25,ytickPOS(iterT),num2str(ytickpos(iterT)),'rotation',0,...
        'horizontalalignment','right','fontsize',14,'interpreter','none')
end

yplot = cat(1,ytickPOS,ytickPOS,NaN(1,numel(ytickPOS)));
xplot = cat(1,zeros(1,numel(ytickPOS))-0.10,zeros(1,numel(ytickPOS))+0.1,zeros(1,numel(ytickPOS)));
plot(xplot,yplot,'k','linewidth',2)
plot([0.1 0.1],[ytickPOS(1) ytickPOS(end)],'k','linewidth',2)


text(xtickpos(1)-(0.5)*range(xtickpos),0,'mm','rotation',90,...
    'color',fontC,'horizontalalignment','right','fontsize',20,'interpreter','none');
pdfName = 'DL_looming_nonJumper_longitudinal_motion.pdf';
text(xtickpos(1),yMax+1,pdfName(1:end-4),...
    'horizontalalignment','left','interpreter','none',...
    'rotation',0,'color',fontC,'fontsize',20);
% {'elevation','aximuth','startSize','stopSize','lv'};
% for iterLoom = 1:4
%     text(xtickpos(end)+1,yMax-iterLoom,pdfName(1:end-4),...
%         'horizontalalignment','left','interpreter','none',...
%         'rotation',0,'color',fontC,'fontsize',20);
% end
if saveFig
    addpath('C:\Users\williamsonw\Documents\MATLAB\ojwoodford-export_fig-165dc92')
    writePath = fullfile('Y:\Ming_RubinLab',pdfName);
    export_fig(writePath,'-nocrop')
    writePath = fullfile('Y:\Ming_RubinLab',[pdfName(1:end-4) '.eps']);
    export_fig(writePath,'-nocrop')
end
