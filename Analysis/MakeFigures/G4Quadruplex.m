%% show that not all late replicating regions are underreplicated in metaphase
load('~/Develop/Mendoza__ReplicationEvolution/Data/DS_stat__200bp.mat');
GC = readtable('~/Develop/Mendoza__ReplicationEvolution/Data/ExternalData/genome_200bpwindows_GC.bed','FileType','text','ReadVariableNames',false);
GC.Properties.VariableNames = {'chr' 'start_point' 'end_point' 'GC'};
T = join(dataset2table(DS) , GC , 'Key' , {'chr' 'start_point' 'end_point'});
T = sortrows(T,'dist_to_the_end_kb','ascend');
%cd('~/Desktop/G4')

% calc DM
%T = T( ~isinf(T.percent_unreplicated_not_trimmed_cdc20) ,:) ; 
DM_window_size = 200;
S = CalcDM_Newman06( T.dist_to_the_end_kb , T.percent_unreplicated_not_trimmed_cdc20_smooth  , DM_window_size , 0 );
T.DM_ypred = S.DM_ypred ; 
T.DM_diff = S.DM_diff ; 

%% load G4 quadruplex data

G4H = readtable( '~/Develop/Mendoza__ReplicationEvolution/Data/ExternalData/Bedrat16/S1_W25_sum_200.bed' ,'Delimiter','\t','FileType','text','Format','%s%d%d%f','TreatAsEmpty' ,'.');
G4H.Properties.VariableNames = {'chr' 'start_point' 'end_point' 'value'};
Q = innerjoin( G4H , T , 'Key', {'chr' 'start_point' 'end_point'});
%%
for I = 1:length(Q)
    if isnan(Q.value(I))
        Q.value(I) = 0;
    end
end
%% Figure
legend_titles = {'low' , 'high'};
clrs1 = lines(6); clrs2 = parula(12); clrs_set = [clrs1(4,:) ; clrs2(6,:)];
fh = figure('units','centimeters','position',[5 5 8 8 ]); 
idx = Q.value > 100 ; 
subplot(1,2,2); hold on; grid on; set(gca , 'FontSize' , 10);
h1 = boxplot(Q.DM_diff , idx ,'symbol','','color' , [.2 .2 .2]);
h = findobj(gca,'Tag','Box');set(h1 , 'LineWidth' , 1.9);
for j=1:length(h)
    patch(get(h(2-j+1),'XData'),get(h(2-j+1),'YData'), clrs_set(2-j+1,:) ,'FaceAlpha',.5 , 'Display' , legend_titles{j});
end
ylim([-15 20]);
ylabel('under-replication DM');
[~,p] = ttest2(Q.DM_diff(Q.value <= 100) , Q.DM_diff(Q.value > 100) )
set(gca,'xtick',[]);
legend('location', 'SouthOutside');
subplot(1,2,1); hold on; grid on; set(gca , 'FontSize' , 10);
h1 = boxplot(Q.percent_unreplicated_not_trimmed_cdc20_smooth , idx ,'symbol','','color' , [.2 .2 .2]);
h = findobj(gca,'Tag','Box');set(h1 , 'LineWidth' , 1.9);
for j=1:length(h)
    patch(get(h(2-j+1),'XData'),get(h(2-j+1),'YData'), clrs_set(2-j+1,:) ,'FaceAlpha',.5 , 'Display' , legend_titles{j});
end
[~,p] = ttest2(Q.percent_unreplicated_not_trimmed_cdc20_smooth(Q.value <= 100) , Q.percent_unreplicated_not_trimmed_cdc20_smooth(Q.value > 100) )
ylim([-15 70]);
ylabel('% unreplicated cells');
set(gca,'xtick',[]);
legend('location', 'SouthOutside');
print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/Figures/Fig2/2g' , '-r300');


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

