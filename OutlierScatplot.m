function [outlierlist] = OutlierScatplot(x1, y1, thresh1, thresh2, absbool, filtbool, xdirection, ydirection, legend1,color1,...
    printboolean)

% DESCRIPTION
% volcano plot of q-value from anova (y1) vs log ratio of average expression levels (x1) 
% with highlighting of genes above or below x and y thresholds

% OUTPUT
% plot + list of highlighted genes

% INPUT
% - x1: expression log ratio 
% - y1: -log of q-value (FDR corrected)
% - thresh1: threshold x-value for outliers
% - thresh2: threshold y-value for outliers
% - absbool: 'abs' for outliers of both tails of x
% - filtbool: '%' for relative outliers (thresh values in %), else absolute
% outliers (thresh values in absolute values of x and y)
% - xdirection: -1 for most extremes values of x (abs or not), +1 for least extremes
% - ydirection: -1 for most extremes values of y, +1 for least extremes
% - legend1: legend name for outliers
% - color1: color for outliers
% - printboolean: boolean for saving fig to pdf or not

xname='Expression log Ratio';
yname='q.anova';

dirName='path/to/dir'

cmap = colormap('parula');
close(gcf);

% remove nans
nonan1 = find(~isnan(y1)&~isinf(y1));
nonan2 = find(~isnan(x1)&~isinf(x1));
if ~isempty(LOI)
    nonanindex= intersect(intersect(nonan1, nonan2),LOI);
else
    nonanindex= intersect(nonan1, nonan2);    
end
y=y1(nonanindex);
x=x1(nonanindex);


% find outliers for x
if strcmp(absbool, 'abs') % both tails of x
    [B, indx] = sortrows(abs(x),xdirection); 
else
    [B, indx] = sortrows(x,xdirection); % one-tailed
end

if strcmp(filtbool, '%') %relative thresholds 
    indseuilx=floor(length(x)*thresh1);
else %absolute thresholds
    if xdirection==1
        if thresh1<0
            indseuilx=min(find(B>-log10(thresh1)));
        else
            indseuilx=min(find(B>log10(thresh1)));
        end
    else
        indseuilx=min(find(B<log10(thresh1)));
    end
end
seuilx=B(indseuilx);

% find outliers for y
[B, indy] = sortrows(y,ydirection);
if thresh2<1 %relative thresholds 
    indseuily=floor(length(y)*thresh2);
else %absolute thresholds
    if ydirection==1
        indseuily=min(find(B>-log(thresh2/100)));
    else
        indseuily=min(find(B<-log(thresh2/100)));
    end
end
seuily=B(indseuily);

% find outliers for both x and y
yellowindex = intersect(indx(1:indseuilx), indy(1:indseuily));

figure1 = figure('Color',[1 1 1]);
set(figure1,'Resize', 'on','PaperOrientation', 'landscape','PaperUnits','points',...
    'PaperSize', [350 350],...
    'PaperPosition', [0 0 350 350], 'InvertHardcopy','off');
hold on
box on
set(gca,'TickLength', [0 0], 'XColor', 'k', 'YColor', 'k');

% plot outliers in color1
h=[];
lgdLst={};
if ~isempty(yellowindex)
    scatter(x(yellowindex), y(yellowindex),200, 'MarkerFaceColor', color1, 'MarkerEdgeColor', color1, 'Marker', '.');
    h(1)=plot(NaN,NaN, 'Marker', '.', 'MarkerSize', 15, 'Color', color1, 'LineStyle', 'none');
    lgdLst{end+1}=[legend1 ' (' num2str(length(yellowindex)) ')'];
end

% plot non outliers with density-colored scatter plot
[nonyellow, ~]=setdiff([1:length(nonanindex)], yellowindex);
scatplot(x(nonyellow),y(nonyellow), 'circles',sqrt((range(x(nonyellow))/10)^2),100, 5, 1, 10);
hold on

xlim([-1.1*max(abs(x)) 1.1*max(abs(x))]);
ylim([0 1.1*max(abs(y))]);

ybornes=exp(-get(gca, 'ylim'));
lb=-floor(log10(ybornes));
lb=10.^-(lb(1):5:lb(2));
rb=[1 0.05 0.01 lb(lb<0.01)];
set(gca, 'YTick', -log(rb), 'YTickLabel', [strread(num2str(rb(rb>=0.01)), '%s')' strread(num2str(rb(rb<0.01),'%10.0e'), '%s')'])

line([-1.1*max(abs(x)) 1.1*max(abs(x))], [seuily seuily],'Color', 'k','LineStyle', '--', 'LineWidth',.5);
line([seuilx seuilx], get(gca, 'ylim'),'Color','k', 'LineStyle', '--', 'LineWidth',.5);
line([-seuilx -seuilx], get(gca, 'ylim'), 'Color','k','LineStyle', '--', 'LineWidth',.5);

colorbar('off');
xlabel(xname, 'Fontsize', 12, 'FontWeight', 'bold');
ylabel(yname, 'Fontsize', 12, 'FontWeight', 'bold');
if ~isempty(h)
    legd=legend(h,lgdLst,'Location','North');
    set(legd,'EdgeColor',[0 0 0],'Color',[1 1 1]);
end
hold off

outlierlist=nonanindex(yellowindex);

figtitle = [xname ' vs ' yname ' per gene_DIPTETRseq ' legend1 ' ' num2str(-xdirection*thresh2*100) '%' ' ' date]

if printboolean
    currF=pwd;
    cd(dirName)
    print(gcf,  '-painters', '-dpdf',[figtitle  '.pdf']);
    cd(currF);
    close(figure1);
end
