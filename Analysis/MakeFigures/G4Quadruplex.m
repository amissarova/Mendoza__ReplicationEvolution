%% show that not all late replicating regions are underreplicated in metaphase
load('~/Develop/Mendoza__ReplicationEvolution/Data/DS_stat__200bp.mat');
GC = readtable('~/Google Drive/CareyLab/ExternalData/Yeast/genome_200bpwindows_GC.bed','FileType','text','ReadVariableNames',false);
GC.Properties.VariableNames = {'chr' 'start_point' 'end_point' 'GC'};
T = join(dataset2table(DS) , GC , 'Key' , {'chr' 'start_point' 'end_point'});
T = sortrows(T,'dist_to_the_end_kb','ascend');
cd('~/Desktop/G4')
%% calc DM
T = T( ~isinf(T.percent_unreplicated_not_trimmed_cdc20) & ~isnan(T.percent_unreplicated_not_trimmed_cdc20) ,:) ; 

DM_window_size = 100 ;
S = CalcDM_Newman06( T.dist_to_the_end_kb , T.percent_unreplicated_not_trimmed_cdc20  , DM_window_size , 0 );
T.DM_ypred = S.DM_ypred ; 
T.DM_diff = S.DM_diff ; 

%%
fn = dir('*.bed');
for I = 1:numel(fn)
    G4H = readtable( fn(I).name ,'Delimiter','\t','FileType','text','Format','%s%d%d%f','TreatAsEmpty' ,'.');
    G4H = G4H( regexpcmp(G4H.Var1,'^chr[IXV]') , :);
    G4H.Properties.VariableNames = {'chr' 'start_point' 'end_point' 'value'};
    G4H = G4H( ~isnan(G4H.value) ,:);
    Q = innerjoin( G4H , T , 'Key', {'chr' 'start_point' 'end_point'});
    Q = sortrows( Q,{'chr'  'arm' 'dist_to_the_end_kb'} ,'descend');
  
    [c,~] = corr( Q.value ,Q.DM_diff ,'rows','complete');
    [c1,~] = corr( Q.value(Q.dist_to_the_end_kb<100) ,Q.DM_diff(Q.dist_to_the_end_kb<100) ,'rows','complete');
    [c2,~] = corr( Q.value(Q.dist_to_the_end_kb<50) ,Q.DM_diff(Q.dist_to_the_end_kb<50) ,'rows','complete');

    fprintf('%s\t%0.02f\t%0.02f\t%0.02f\n' , fn(I).name , c, c1, c2);
  
%   % I want to know if having lots of G4s along the end of an arm
%   %  is one of the reasons for arm-to-arm differences in max %unrep
%   %  I think this can be ansswered with cumsum() of G4 score along chr arm
%   % but needs more thought. 
%   Q = Q(Q.dist_to_the_end_kb<20,:);
%     Q.ID = strcat( Q.chr , Q.arm);
%     uchr = unique(Q.ID);
%     Q.cumsum = NaN( height(Q),1);
%     for uchrI = 1:numel(uchr)
%         idx = strcmp(Q.ID,uchr{uchrI});
%         Q.cumsum(idx) = cumsum(Q.value(idx));
%     end

%     R = Q(Q.dist_to_the_end_kb > 0 & Q.dist_to_the_end_kb <= 15,:);
%     G = grpstats( R , {'chr' 'arm'} , 'mean' , 'DataVars' , {'percent_unreplicated_not_trimmed_cdc20' 'value' 'cumsum'});
%     
%     figure; 
%     plot( G.mean_percent_unreplicated_not_trimmed_cdc20 , G.mean_cumsum , 'ok')
%     [c,p]=corr(G.mean_percent_unreplicated_not_trimmed_cdc20 , G.mean_cumsum ,'rows','complete')
  % 

end
disp('\n')

%% Figure

G4H = readtable( 'S1_W25_sum_200.bed' ,'Delimiter','\t','FileType','text','Format','%s%d%d%f','TreatAsEmpty' ,'.');
G4H.Properties.VariableNames = {'chr' 'start_point' 'end_point' 'value'};
Q = innerjoin( G4H , T , 'Key', {'chr' 'start_point' 'end_point'});

fh = figure('units','centimeters','position',[5 5 4 6 ]);
idx = Q.value > 100 ; 
boxplot(Q.DM_diff , idx ,'symbol','')
ylim([-15 50])
ylabel('% underreplication')
set(gca,'xtick',[])