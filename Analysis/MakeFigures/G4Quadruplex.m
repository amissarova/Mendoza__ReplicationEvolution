%% show that not all late replicating regions are underreplicated in metaphase
load('~/Develop/Mendoza__ReplicationEvolution/Data/DS_stat__200bp.mat');
GC = readtable('~/Google Drive/CareyLab/ExternalData/Yeast/genome_200bpwindows_GC.bed','FileType','text','ReadVariableNames',false);
GC.Properties.VariableNames = {'chr' 'start_point' 'end_point' 'GC'};
T = join(dataset2table(DS) , GC , 'Key' , {'chr' 'start_point' 'end_point'});
T = sortrows(T,'dist_to_the_end_kb','ascend');
cd('~/Desktop/G4')

%% load G4 quadruplex data
G4H = readtable( '~/Develop/Mendoza__ReplicationEvolution/Data/ExternalData/Bedrat16/S1_W25_sum_200.bed' ,'Delimiter','\t','FileType','text','Format','%s%d%d%f','TreatAsEmpty' ,'.');
G4H.Properties.VariableNames = {'chr' 'start_point' 'end_point' 'value'};
Q = innerjoin( G4H , T , 'Key', {'chr' 'start_point' 'end_point'});
Q.value( isnan(Q.value)) = 0 ; % '.' in bedtools map is from data w/no value. here that's the same as zero
%% Figure
fh = figure('units','centimeters','position',[5 5 4 6 ]);
idx = Q.value > prctile(Q.value,99.99) ; 
boxplot(Q.percent_underreplicated_cdc20_not_trimmed_DM_dist , idx ,'symbol','')
ylim([-15 50])
ylabel('% underreplication')
set(gca,'xtick',[])


%% playing with more analysis down here
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

