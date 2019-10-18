%% Is the segregation of rDNA and telomeric loci sensitive to HU? 
% reviewer question: 
%   2. Also related to DAPI-poor areas, if the nucleolus is a major site for late mitotic replication, does polymerase inhibition in mitosis preferentially affect rDNA segregation?
% load data
DD = '~/Develop/Mendoza__ReplicationEvolution/Data/' ; 
rDNA  = readtable( [DD 'rdna segregation from anaphase -- October 2019.xlsx'] ,'Sheet','rDNA');
tel12 = readtable( [DD 'rdna segregation from anaphase -- October 2019.xlsx'] ,'Sheet','tel12');
tel12S = stack(tel12,tel12.Properties.VariableNames) ; 
rDNAS = stack(rDNA  , rDNA.Properties.VariableNames) ;

%% histograms

fh = figure('units','centimeters','position',[5 5 8 8]);
hold on ;
[t,c]=count_unique(rDNA.untreated);
xl = 0.9:4 ; 
bar(xl,c,0.4,'FaceColor',[.4 .4 .4])
for I = 1:numel(t)
    v = bootstrp( 10000 , @(X)sum(X==t(I)) , rDNA.untreated)  ; 
    errorbar( xl(I) , mean(v) , std(v),'.k')
end

[t,c]=count_unique(rDNA.HU);
xl = 1.2:5;
bar(xl,c,0.4,'FaceColor',clrs(2,:))
for I = 1:numel(t)
    v = bootstrp( 10000 , @(X)sum(X==t(I)) , rDNA.HU)  ; 
    errorbar( xl(I) , mean(v) , std(v),'.','Color',clrs(2,:))
end

set(gca,'xtick',1:4)
set(gca,'xticklabel',{'0 min' '4 min' '8 min' '12 min'})
xlim([0.6 4.5])
ylim([0 115])
ylabel('# of cells')
xlabel('Time to segregation from anaphase onset')
[~,p]=ttest2(rDNA.HU,rDNA.untreated)

%% tel
fh = figure('units','centimeters','position',[5 5 8 8]);
hold on ;
[t,c]=count_unique(tel12.untreated);
xl = 0.9:5 ; 
bar(xl,c,0.4,'FaceColor',[.4 .4 .4])
for I = 1:numel(t)
    v = bootstrp( 10000 , @(X)sum(X==t(I)) , tel12.untreated)  ; 
    errorbar( xl(I) , mean(v) , std(v),'.k')
end

[t,c]=count_unique(tel12.HU); t = vertcat(0,t); c = vertcat(0,c) ; 
xl = 1.2:6;
bar(xl,c,0.4,'FaceColor',clrs(2,:))
for I = 1:numel(t)
    v = bootstrp( 10000 , @(X)sum(X==t(I)) , tel12.HU)  ; 
    errorbar( xl(I) , mean(v) , std(v),'.','Color',clrs(2,:))
end

set(gca,'xtick',1:5)
set(gca,'xticklabel',{'0 min' '4 min' '8 min' '12 min' '16 min'})
xlim([0.6 5.5])
ylim([0 94])
ylabel('# of cells')
xlabel('Time to segregation from anaphase onset')
[~,p]=ttest2(tel12.HU,tel12.untreated)


%% boxplots show no difference in the distribution

figure;
hold on ;
histogram(tel12.untreated)
histogram(tel12.HU)
legend({'untreated' 'HU'})
title('subtel 12R')
xlabel('Minutes to segregation from anaphase onset')
ylabel('# of cells')
[~,p]=ttest2(tel12.HU,tel12.untreated)

x=6;
x = [sum(tel12.HU<x) sum(tel12.HU>x)  ; sum(tel12.untreated<x) sum(tel12.untreated>x)  ] ;
[~,p,stats]=fishertest(x)


% Main figure; 
figure('units','centimeters','position',[5 5 4 8]); 
clrs = get(gca,'ColorOrder');
clrs = [0 0 0 ; clrs(2,:)] ;
bh = boxplot(tel12S.HU_untreated , tel12S.HU_untreated_Indicator,'symbol','','Color','k','GroupOrder',{'untreated' 'HU'});
set(bh,'LineWidth',1.3);
h = findobj(gca,'Tag','Box');
ylim([0 15])
K = numel(h);
for j = 1:K
        patch(get(h(K-j+1),'XData'),get(h(K-j+1),'YData'), clrs(j,:) ,'FaceAlpha',.4 );
end    
ylabel('Minutes to segregation from anaphase onset')
%set(gca,'xticklabel',[])

[~,p]=ttest2(tel12.HU,tel12.untreated)
title('subtel 12R')

figure('units','centimeters','position',[5 5 4 8]); 
clrs = get(gca,'ColorOrder');
clrs = [0 0 0 ; clrs(2,:)] ;
bh = boxplot(rDNAS.HU_untreated , rDNAS.HU_untreated_Indicator,'symbol','','Color','k','GroupOrder',{'untreated' 'HU'});
set(bh,'LineWidth',1.3);
h = findobj(gca,'Tag','Box');
ylim([0 15])
K = numel(h);
for j = 1:K
        patch(get(h(K-j+1),'XData'),get(h(K-j+1),'YData'), clrs(j,:) ,'FaceAlpha',.4 );
end    
ylabel('Minutes to segregation from anaphase onset')
%set(gca,'xticklabel',[])

[~,p]=ttest2(rDNA.HU,rDNA.untreated)
title('rDNA')