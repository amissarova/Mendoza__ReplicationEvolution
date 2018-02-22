%% load data
figname = 'DMcorrected_vs_UnCorr___Trep_Dist2ARS' ; 
FigHeight = 8 ; 
load('~/Develop/Mendoza__ReplicationEvolution/Data/DS_stat__200bp.mat');
GC = readtable('~/Google Drive/CareyLab/ExternalData/Yeast/genome_200bpwindows_GC.bed','FileType','text','ReadVariableNames',false);
GC.Properties.VariableNames = {'chr' 'start_point' 'end_point' 'GC'};
T = join(dataset2table(DS) , GC , 'Key' , {'chr' 'start_point' 'end_point'});
idx_keep = ~isnan(T.percent_unreplicated_not_trimmed_cdc20) & ~isnan(T.percent_unreplicated_not_trimmed_dbf2) ...
    & ~isinf(T.percent_unreplicated_not_trimmed_cdc20) & ~isinf(T.percent_unreplicated_not_trimmed_dbf2) ; 

T = T( idx_keep ,:);
T = sortrows( T , 'dist_to_the_end_kb' ,'ascend');
%% calculate DM
DM_window_size = 1e2 ; % size in # of 200bp windows
Scdc20 = CalcDM_Newman06( T.dist_to_the_end_kb , T.percent_unreplicated_not_trimmed_cdc20  , DM_window_size , 0 );
Sdbf2 = CalcDM_Newman06( T.dist_to_the_end_kb , T.percent_unreplicated_not_trimmed_dbf2 ,  DM_window_size , 0 );


%% plot DM Cartoon
fh = figure('units','centimeters','position',[5 5 16 FigHeight ]); 
hold on ; 
sh = scatter( Scdc20.X , Scdc20.Y , 20 , 'k' ,'MarkerFaceColor',[.7 .7 .7],'MarkerEdgeAlpha',0.15,'MarkerFaceAlpha',0.15);
plot( Scdc20.X , Scdc20.DM_ypred ,'-r')
xlim([0.1 500])
ylim([-19 75])
grid on ;
xlabel('KB to the chromosome end')
ylabel('% underreplicated in metaphase')
print('-dpng',[figname '_linX.png'],'-r300')

xlim([0 105])
print('-dpng',[figname '_linZoom.png'],'-r300')

xlim([0.1 500])
set(gca,'xscale','log')
set(gca,'xtick',[0.2 0.5 1 2 5 10 25 100 500])
print('-dpng',[figname '_logX.png'],'-r300')

close ; 

%%
fh = figure('units','centimeters','position',[5 5 8 FigHeight ]); 
hold on ; 
sh = scatter( Scdc20.X , Scdc20.DM_diff , 20 , 'k' ,'MarkerFaceColor',[.7 .7 .7],'MarkerEdgeAlpha',0.15,'MarkerFaceAlpha',0.15);
line(xlim,[0 0],'LineStyle','-','Color','r','LineWidth',2)
xlim([0 500])
ylim([-49 49])
grid on ;
xlabel('KB to the chromosome end')
ylabel('Underreplication DM (residuals)')
set(gca,'xscale','lin')
%set(gca,'xtick',[0.2 0.5 1 2 5 10 25 100 500])
print('-dpng',[figname '_DMresid.png'],'-r300')
close ; 

%% plot % underrep vs Dist 2 ARS
fh = figure('units','centimeters','position',[5 5 8 FigHeight ]); 
hold on ; 

clrs = get(gca,'ColorOrder');
clr_uncorr = [.7 .7 .7] ; 
clr_corr   = [0 0 0]  ; 

xl = round(T.dist_to_ARS./10)*10 ; 
uxl = unique(xl);

bh1 = boxplot( T.percent_unreplicated_not_trimmed_cdc20 , xl , 'plotstyle','compact','symbol','','Positions',uxl-0.5,'Color',clr_uncorr,'labelorientation','horizontal') ;
%bh2 = boxplot( Scdc20.Exp_1_diff , xl , 'plotstyle','compact','symbol','','Positions',uxl+1,'Color',clrs(1,:),'labelorientation','horizontal')
%bh2 = boxplot( Scdc20.Exp_2_diff , xl , 'plotstyle','compact','symbol','','Positions',uxl+2,'Color',clrs(2,:),'labelorientation','horizontal')
bh3 = boxplot( Scdc20.DM_diff , xl , 'plotstyle','compact','symbol','','Positions',uxl+1,'Color',clr_corr,'labelorientation','horizontal') ;
ylim([-5 15])
%arrayfun( @(X) set(X,'MarkerFaceColor',clr_uncorr) , bh1(:))
%arrayfun( @(X) set(X,'MarkerFaceColor',clr_corr) , bh3(:))
xlabel('KB to the nearest ARS')
ylabel('% underreplicated in metaphase')
fprintf('DM: %0.02f\nEx1 %0.02f\nEx2 %0.02f' , mean(abs(Scdc20.DM_diff )), mean(abs(Scdc20.Exp_1_diff )) ,  mean(abs(Scdc20.Exp_2_diff )));
[~,p] = ttest2( Scdc20.DM_diff(xl==35) , Scdc20.DM_diff(xl~=35))


%% plot % underrep vs Trep
fh = figure('units','centimeters','position',[5 5 9 FigHeight ]); 
hold on ; 

clrs = get(gca,'ColorOrder');
clr_uncorr = [.7 .7 .7] ; 
clr_corr   = [0 0 0]  ; 

xl = round(T.Trep_spline./10)*10 ;
xl(xl<10)=10;
%xl(xl>60)=60 ;
uxl = unique(xl);

bh1 = boxplot( T.percent_unreplicated_not_trimmed_cdc20 , xl , 'plotstyle','compact','symbol','','Positions',uxl-0.5,'Color',clr_uncorr,'labelorientation','horizontal') ;
%bh2 = boxplot( Scdc20.Exp_1_diff , xl , 'plotstyle','compact','symbol','','Positions',uxl+1,'Color',clrs(1,:),'labelorientation','horizontal')
%bh2 = boxplot( Scdc20.Exp_2_diff , xl , 'plotstyle','compact','symbol','','Positions',uxl+2,'Color',clrs(2,:),'labelorientation','horizontal')
bh3 = boxplot( Scdc20.DM_diff , xl , 'plotstyle','compact','symbol','','Positions',uxl+1,'Color',clr_corr,'labelorientation','horizontal') ;
ylim([-5 15])
line(xlim,[0 0],'Color',[.5 .5 .5])

%arrayfun( @(X) set(X,'MarkerFaceColor',clr_uncorr) , bh1(:))
%arrayfun( @(X) set(X,'MarkerFaceColor',clr_corr) , bh3(:))
xlabel('Replication timing (min after \alpha factor release)')
ylabel('% underreplicated in metaphase')
fprintf('DM: %0.02f\nEx1 %0.02f\nEx2 %0.02f' , mean(abs(Scdc20.DM_diff )), mean(abs(Scdc20.Exp_1_diff )) ,  mean(abs(Scdc20.Exp_2_diff )));

[ uxl arrayfun( @(X) ttest2( Scdc20.DM_diff(xl==X) , Scdc20.DM_diff ,'tail','right' , 'alpha' , 0.05./numel(uxl)) , uxl) ]

%% plot %GC underrep vs Trep
fh = figure('units','centimeters','position',[5 5 9 FigHeight ]); 
hold on ; 

clrs = get(gca,'ColorOrder');
clr_uncorr = [.7 .7 .7] ; 
clr_corr   = [0 0 0]  ; 

xl = round(T.GC*20)./20 ; 
%xl(xl>60)=60 ;
uxl = unique(xl);

bh1 = boxplot( T.percent_unreplicated_not_trimmed_cdc20 , xl , 'plotstyle','compact','symbol','','Positions',uxl-0.005,'Color',clr_uncorr,'labelorientation','horizontal') ;
%bh2 = boxplot( Scdc20.Exp_1_diff , xl , 'plotstyle','compact','symbol','','Positions',uxl+1,'Color',clrs(1,:),'labelorientation','horizontal')
%bh2 = boxplot( Scdc20.Exp_2_diff , xl , 'plotstyle','compact','symbol','','Positions',uxl+2,'Color',clrs(2,:),'labelorientation','horizontal')
bh3 = boxplot( Scdc20.DM_diff , xl , 'plotstyle','compact','symbol','','Positions',uxl+0.005,'Color',clr_corr,'labelorientation','horizontal') ;
ylim([-9 60])
%arrayfun( @(X) set(X,'MarkerFaceColor',clr_uncorr) , bh1(:))
%arrayfun( @(X) set(X,'MarkerFaceColor',clr_corr) , bh3(:))
xlabel('Percent GC')
ylabel('% underreplicated in metaphase')
fprintf('DM: %0.02f\nEx1 %0.02f\nEx2 %0.02f' , mean(abs(Scdc20.DM_diff )), mean(abs(Scdc20.Exp_1_diff )) ,  mean(abs(Scdc20.Exp_2_diff )));
xt = get(gca,'xtick');
xtl = get(gca,'xticklabel');

set(gca,'xtick',xt(1:2:end))
set(gca,'xticklabel',xtl(1:2:end))

sigmat = [ uxl arrayfun( @(X) ttest2( Scdc20.DM_diff(xl==X) , Scdc20.DM_diff ,'tail','right' , 'alpha' , 0.05./numel(uxl)) , uxl) ]
text( uxl(find(sigmat(:,2))) , repmat(55,sum(sigmat(:,2)),1) ,'*')
line(xlim,[0 0],'Color',[.5 .5 .5])
