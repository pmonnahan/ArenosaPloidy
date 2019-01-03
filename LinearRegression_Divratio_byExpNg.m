function [mdlTot] = LinearRegression_Divratio_byExpNg(RNAseqNDE, foldPos, POPtetStats, POPdipStats, Ng, yname,...
    threshd,threshu,threshd2,threshu2,expthreshd, expthreshu, pithreshu, TetArray, DipArray,diptet_subset, AlTranscriptIDnum, HTcounts, filename, printBool )

% DESCRIPTION
% multiple linear model of 0/4 ratio of theta estimates of diversity by ploidy,
% expression, and effective Ng

% INPUT
% - RNAseqNDE: set of genes defined as NDE (non-differentially expressed)
% between 2x and 4x
% - foldPos: the columns index of all/0-dg/4-dg theta estimates by gene within
% geneStats matrix (the respective nb of sites for each estimate are in
% the same order three columns to the left in the geneStats matrix)
% - POPtetStats: geneStats matrix for tetraploid populations
% - POPdipStats: geneStats matrix for diploid populations
% - Ng: estimates of effective nbs of haploid genomes per population
% - yname:  name of estimator of theta used ('theta_W' or 'theta_pi')
% - threshd,threshu: lower and upper threshold for nb of 0-dg sites required by gene 
% - threshd2,threshu2: lower and upper threshold for nb of 4-dg sites required by gene
% - expthreshd, expthreshu: lower and upper threshold of linear range of
% theta by expression
% - pithreshu: upper threshold for diversity
% - TetArray: cell array of tetraploid populations names
% - DipArray: cell array of diploid populations names
% - diptet_subset: name of subset of populations to run analysis on (used to
% remove admixed tetraploid populations or restrict the analysis to w-carp
% populations) ('ALL', 'no_admX', or 'W_Carp_only')
% - AlTranscriptIDnum: A. lyrata gene IDs
% - HTcounts: Structure containing gene expression data in the order of AlTranscriptIDnum
% - filename: name to give to figure
% - printBool: boolean to print figure to pdf or not


dirName='path/to/dir'

close gcf
figure1 = figure('Color',[1 1 1]);
cmap=colormap('parula');
set(figure1,'Resize', 'on','PaperOrientation', 'landscape','PaperUnits',...
    'points','PaperSize', [1050 600], 'PaperPosition', [0 0 1050 600],...
    'InvertHardcopy','off');
hold on
samplename='';
legendIND=[1 3; 2 4];
lineType={'-' '--'};
foldName={'' '0fold' '4fold'};
legendNames=cell(6,1);
Nploid=5;
idip=[1 2 3 4 5 6 7 8 9];
Ndip=length(idip);
itet=[1 2 3 4 5 6 7 8 9];
Ntet=length(itet);

colorPop={[.7 .8 .9]  [.7 .9 .8]};
colorAll= {cmap(5,:) [.7 .8 .9] ;cmap(35,:)  [.7 .9 .8]};

PopOI={'SNO' 'SCH'}; % outliers
WCPopOI={'HNI' 'SNO' 'TRD' 'VEL' 'SPI' 'TKO' 'TRT' 'ZAP' 'TRE'}; %w-carp. pops
AdmXPopOI={'STE' 'TBG' 'KOW' 'DRA' 'LAC' 'TZI'}; %admixed pops
xrate=0.5;

subplot(2,2,1)
hold on

Xmin = 0;
Xmax = 15000;
expressionVecTot=[]; diversityVecTot={[] [] []}; NgVecTot=[]; ploidyVecTot={}; gNameVecTot=[]; popNameVecTot=[]; 

for iploid=1:Ndip+Ntet
    if iploid<=Ndip
        reIND=iploid;
        maxIND=Ndip;
        popStats = POPdipStats{idip(iploid)}.geneStats;
        popname = DipArray{idip(iploid)};
        ploidyVec ={'diploid'};
        ploidIND=1; 
        RNAseq.Genes=popStats(:,1);
        [~, i1, i2] = intersect(RNAseq.Genes, AlTranscriptIDnum);
        RNAseq.expression=HTcounts.Expression.Dip.Means(i2);
    else
        reIND=iploid-Ndip;
        maxIND=Ntet;
        popStats = POPtetStats{itet(reIND)}.geneStats;
        popname = TetArray{itet(reIND)};
        ploidyVec ={'tetraploid'};
        ploidIND=2;
        RNAseq.Genes=popStats(:,1);
        [~, i1, i2] = intersect(RNAseq.Genes, AlTranscriptIDnum);
        RNAseq.Genes=RNAseq.Genes(i1);
        RNAseq.expression=HTcounts.Expression.Tet.Means(i2);
    end
    
    % filter out genes with nb of 0-dg sites outside of [threshd - threshu] and
    % 4-dg sites outside of [threshd2 - threshu2] and 0-dg theta > pithreshu
    % and null expression levels
    NonNULLval1 = find(((popStats(i1,foldPos(2))==0)+(popStats(i1,foldPos(2)-3)<threshd)+(popStats(i1,foldPos(2)-3)>threshu)+...
        (popStats(i1,foldPos(3))==0)+(popStats(i1,foldPos(3)-3)<threshd2)+(popStats(i1,foldPos(3)-3)>threshu2)...
        +(popStats(i1,foldPos(2))>pithreshu)+(RNAseq.expression==0))==0); 

    [~, IntVal, ~]= intersect(popStats(i1,1),RNAseqNDE.genes); %remove genes not NDE (non-differentially expressed between 2x vs 4x)
    NonNULLval=intersect(NonNULLval1,IntVal);
    LTet=length(NonNULLval);
    genenameVec=(RNAseq.Genes(NonNULLval(:)));
    expressionVec=(RNAseq.expression(NonNULLval(:)));
    diversityVec=log(popStats(i1(NonNULLval(:)),foldPos(2))./popStats(i1(NonNULLval(:)),foldPos(3)));
    diversitySVec=(popStats(i1(NonNULLval(:)),foldPos(3)));
    Xmin=max(Xmin, min(expressionVec));
    Xmax=min(Xmax, max(expressionVec));
    
    if strcmp('ALL', diptet_subset)
        gNameVecTot=vertcat(gNameVecTot, genenameVec);
        NgVecTot=vertcat(NgVecTot, Ng(iploid)*ones(LTet,1));
        expressionVecTot=vertcat(expressionVecTot, expressionVec);
        diversityVecTot{2}=vertcat(diversityVecTot{2}, diversityVec);
        diversityVecTot{3}=vertcat(diversityVecTot{3}, diversitySVec);
        ploidyVecTot=vertcat(ploidyVecTot, ploidyVec(ones(LTet,1)));
        popnameVec ={popname};
        popNameVecTot=vertcat(popNameVecTot, popnameVec(ones(LTet,1)));

    elseif strcmp('W_Carp_only', diptet_subset)
        if max(strcmp(WCPopOI, popname)) % RUN MODEL ON West CARP only
            popname
            gNameVecTot=vertcat(gNameVecTot, genenameVec);
            NgVecTot=vertcat(NgVecTot, Ng(iploid)*ones(LTet,1));
            expressionVecTot=vertcat(expressionVecTot, expressionVec);
            diversityVecTot{2}=vertcat(diversityVecTot{2}, diversityVec);
            diversityVecTot{3}=vertcat(diversityVecTot{3}, diversitySVec);
            ploidyVecTot=vertcat(ploidyVecTot, ploidyVec(ones(LTet,1)));
            popnameVec ={popname};
            popNameVecTot=vertcat(popNameVecTot, popnameVec(ones(LTet,1)));
        end
        
    elseif strcmp('no_admX', diptet_subset)
        if ~max(strcmp(AdmXPopOI, popname)) % Exclude ADMIXED TETs
            popname
            gNameVecTot=vertcat(gNameVecTot, genenameVec);
            NgVecTot=vertcat(NgVecTot, Ng(iploid)*ones(LTet,1));
            expressionVecTot=vertcat(expressionVecTot, expressionVec);
            diversityVecTot{2}=vertcat(diversityVecTot{2}, diversityVec);
            diversityVecTot{3}=vertcat(diversityVecTot{3}, diversitySVec);
            ploidyVecTot=vertcat(ploidyVecTot, ploidyVec(ones(LTet,1)));
            popnameVec ={popname};
            popNameVecTot=vertcat(popNameVecTot, popnameVec(ones(LTet,1)));
        end
    end
end


% 3 EXPRESSION ranges : < expthreshd ; > expthreshd & < expthreshu ; > expthreshu
rangeLims={[Xmin expthreshd]; [expthreshd expthreshu]; [expthreshu Xmax]};
for reg=2 %run regression within linear range (within [expthreshd expthreshu])
    rangeName{reg}=[num2str(rangeLims{reg}(1)) '-' num2str(rangeLims{reg}(2))];
    indofRange=find((expressionVecTot<rangeLims{reg}(1))+(expressionVecTot>rangeLims{reg}(2))==0);
    lengthofRange{reg}=length(indofRange);
    
    %standardize predictors to compare effect estimates    
    log_Exp_vec_std= zscore(log(expressionVecTot(indofRange)));
    Ng_vec_std= zscore(NgVecTot(indofRange));  
    
    %with standardized predictors 
    dsTot=dataset( ploidyVecTot(indofRange),log_Exp_vec_std, diversityVecTot{2}(indofRange), Ng_vec_std,...
        'Varnames',{'ploidy', 'logExpression','theta_ratio', ' Ng'});
    
    % export data to txt file
    currdir=pwd;
    cd('C:\Users\pbaduel\Dropbox\JIC\Labing\Tables')
    export(dsTot,'File',['Sample ' filename ' ' yname ' vs ' 'Expression' ' resp_ploidy_thS_Ng '...
    'NDE' 'allpop ' num2str(threshd) ' SNPfilter exp' num2str(expthreshd) '-' num2str(expthreshu) 'maxfilter ' date '.txt'])
    cd(currdir)
    
    %multiple linear model with robust estimators
    mdlTot = fitlm(dsTot,'theta_ratio ~ 1 + logExpression + ploidy  + Ng + logExpression*ploidy ','RobustOpts','on')
     
    %store text output of linear regression
    mdlOutput{reg}=evalc('disp(mdlTot)');
    mdlOutput{reg} = strrep(mdlOutput{reg},'<strong>','\bf');
    mdlOutput{reg} = strrep(mdlOutput{reg},'</strong>','\rm');
    mdlOutput{reg} = strrep(mdlOutput{reg},'_','\_');
   
    %plot estimated effects
    h = plotEffects(mdlTot)
    h(1).Marker='+';
    h(1).MarkerSize=2;
    h(1).LineWidth=1;
    h(1).MarkerEdgeColor=[0 0 0];
    for j=2:4
        h(j).Color=[0 0 0];
        h(j).LineWidth=1;
    end
end
    
    
posVec1=get(gca,'position')
set(gca, 'TickLength',[0 0], 'ycolor','k', 'xcolor', 'k'); 
box on
hold off

%plot interaction terms between Ng and logExpression
subplot(2,2,3)
hold on
h=plotInteraction(mdlTot,'Ng','logExpression')
for j=[2 3]
    h(j).Color=[0 0 0];
    h(j).LineWidth=1;
end
for j=[4 8]
    h(j).Marker='+';
    h(j).LineWidth=1;
    h(j).MarkerSize=2;
    h(j).MarkerEdgeColor=[0.5 0.5 0.5];
end
for j=[5 6 7 9 10 11]
    h(j).Color=[0.5 0.5 0.5];
    h(j).LineWidth=1;
end
h(1).Marker='+';
h(1).LineWidth=1;
h(1).MarkerSize=2;
h(1).MarkerEdgeColor=[0 0 0];
 
set(gca, 'TickLength',[0 0], 'ycolor','k', 'xcolor', 'k'); 
posVec2=get(gca,'position')
box on
hold off

%plot interaction terms between logExpression and ploidy
subplot(2,2,4)
hold on
h=plotInteraction(mdlTot,'logExpression','ploidy')
for j=[2 3]
    h(j).Color=[0 0 0];
    h(j).LineWidth=1;
end
for j=[4 7]
    h(j).Marker='+';
    h(j).LineWidth=1;
    h(j).MarkerSize=2;
    h(j).MarkerEdgeColor=[0.5 0.5 0.5];
end
for j=[5 6 7 8 9 10]
    h(j).Color=[0.5 0.5 0.5];
    h(j).LineWidth=1;
end
h(1).Marker='+';
h(1).LineWidth=1;
h(1).MarkerSize=2;
h(1).MarkerEdgeColor=[0 0 0];
 
set(gca, 'TickLength',[0 0],'xlim',[-0.25 .1], 'ycolor','k', 'xcolor', 'k'); 
posVec2=get(gca,'position')
box on
hold off 
 
%print text output of linear model
for reg=2
    subplot(2,2,2)
    hold on
    xlim([0 1])
    ylim([0 1])
    axis off    

    FixedWidth = get(0,'FixedWidthFontName');
    text(-0.1, 1, [rangeName{reg} ' ' mdlOutput],'FontName',FixedWidth,...
        'Fontsize', 8,'Interpreter','Tex', 'HorizontalAlignment', 'left','VerticalAlignment', 'top');
     hold off
end

hold all
title(diptet_subset, 'interpreter', 'none')

%save fig
figtitle = ['Sample ' filename '-' diptet_subset ' ' ['0-4 ratio  ' yname] ' vs ' 'Expression' ' resp_ploidy_thS_Ng '...
    'NDE' ' model effects ' downsampling 'allpop ' num2str(threshd) ' SNPfilter exp' num2str(expthreshd) '-' num2str(expthreshu) 'maxfilter ' date]

if printBool
    currF=pwd;
    cd(dirName)
    print(gcf,  '-painters', '-dpdf',[figtitle '.pdf']);
    close(gcf);
    cd(currF);
end
