%% 2019.10.17 -- new data on nuclear segregation and chr bridges
% Figure1_nuclear_segragation_and_bridges_new_data.m


%% load data
PROJECTDIR = '~/CareyLab/Projects/2017__MendozaReplication/' ;
DATADIR = '~/Develop/Mendoza__ReplicationEvolution/2019.10.17 -- new data on nuclear segregation and chr bridges/' ; 
DATA_FILE_1 = [ DATADIR 'psf2 bridge lifetime pooled.txt' ] ;
DATA_FILE_2 = [ DATADIR 'HU bridge lifetime pooled.txt' ] ;

FIGUREOUTDIR = [ PROJECTDIR '/Figures_for_revision/' ] ;  system(['mkdir -p ' FIGUREOUTDIR]);

T1 = readtable( DATA_FILE_1 ) ; 
T1.Properties.VariableNames = {'WT' 'psf2' 'rad9' 'psf2rad9'} ; 

T2 = readtable( DATA_FILE_2 ) ; 
T2.Properties.VariableNames = {'WT' 'HU' 'rad9' 'rad9HU'} ; 

S1 = stack(T1,1:4) ; 
S2 = stack(T2,1:4) ; 
S1.Properties.VariableNames = {'condition' 'minutes'};
S2.Properties.VariableNames = {'condition' 'minutes'};

%% define colors for the lines
figure; 
clrs = [0 0 0 ; get(gca,'ColorOrder') ] ;
clrs1 = clrs( [1:2 4:end] , :);
clrs2 = [ clrs( [1 3 ] , :) ; [.5 .5 .5] ; clrs(4,:) ] ;
clrs1 = [ clrs( [1 2 ] , :) ; [.5 .5 .5] ; clrs(5,:) ] ;
close ; 
%% Main figure; 
figure('units','centimeters','position',[5 5 6 5]); 

bh = boxplot(S1.minutes , S1.condition,'symbol','','notch','on','Color',[.2 .2 .2]);
set(bh,'LineWidth',1.3);
h = findobj(gca,'Tag','Box');
ylim([0 20])

K = numel(h);
for j = 1:K
        patch(get(h(K-j+1),'XData'),get(h(K-j+1),'YData'), clrs1(j,:) ,'FaceAlpha',.8 );
end    
ylabel('Minutes after release')
set(gca,'xticklabel',[])

[~,p]=ttest2(T1.WT,T1.psf2)

%% Main figure; 
figure('units','centimeters','position',[5 5 6 5]); 

bh = boxplot(S2.minutes , S2.condition,'symbol','','notch','on','Color',[.2 .2 .2]);
set(bh,'LineWidth',1.3);
h = findobj(gca,'Tag','Box');
ylim([0 20])

K = numel(h);
for j = 1:K
        patch(get(h(K-j+1),'XData'),get(h(K-j+1),'YData'), clrs2(j,:) ,'FaceAlpha',.8 );
end    
ylabel('Minutes after release')
set(gca,'xticklabel',[])

[~,p]=ttest2(T2.WT,T2.HU)
