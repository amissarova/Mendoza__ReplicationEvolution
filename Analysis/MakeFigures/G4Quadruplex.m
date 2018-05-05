%% show that not all late replicating regions are underreplicated in metaphase
load('~/Develop/Mendoza__ReplicationEvolution/Data/DS_stat__200bp.mat');
GC = readtable('~/Google Drive/CareyLab/ExternalData/Yeast/genome_200bpwindows_GC.bed','FileType','text','ReadVariableNames',false);
GC.Properties.VariableNames = {'chr' 'start_point' 'end_point' 'GC'};
T = join(dataset2table(DS) , GC , 'Key' , {'chr' 'start_point' 'end_point'});
T = sortrows(T,'dist_to_the_end_kb','ascend');


%% load G4 quadruplex data from my predictions
G4H = readtable( '~/Develop/Mendoza__ReplicationEvolution/Data/ExternalData/Bedrat16/S1_W25_sum_200.bed' ,'Delimiter','\t','FileType','text','Format','%s%d%d%f','TreatAsEmpty' ,'.');
G4H.Properties.VariableNames = {'chr' 'start_point' 'end_point' 'value'};
Q = innerjoin( G4H , T , 'Key', {'chr' 'start_point' 'end_point'});
Q.value( isnan(Q.value)) = 0 ; % '.' in bedtools map is from data w/no value. here that's the same as zero
Q = sortrows(Q,'value','descend');

%% load G4 quad data from Capra 2010
Capra10 = readtable('~/Data/Capra10/Dataset_S1.TXT');
Capra10.Var2 = [];
Capra10.Var3 = [];
Capra10.Var6 = [];
Capra10.Var8 = [];
Capra10.Var9 = [];
Capra10.Var1 = str2double(regexprep( Capra10.Var1 , 'scer_chr',''));
Capra10.Properties.VariableNames = {'chr_num' 'start_point' 'end_point' 'strand'} ;
Capra10.CapraPos = round((Capra10.end_point - Capra10.start_point )/2) +   Capra10.start_point ;

T.mid_point = T.start_point + 100 ;
T.Capra10_dist = NaN( height(T)  , 1);
for I = 1:height(T)
    idx = Capra10.chr_num == T.chr_num(I);
    T.Capra10_dist(I) = min(abs(Capra10.start_point(idx) - T.mid_point(I))) ; 
end
%% %% Capra10 boxplot figure
idxnan = isnan(T.percent_underreplicated_cdc20_not_trimmed_DM_dist) | isinf(T.percent_underreplicated_cdc20_not_trimmed_DM_dist) ;
fh = figure('units','centimeters','position',[5 5 20 20 ]);
hold on ;
idx =  T.Capra10_dist < 50 ;
histogram( T.percent_underreplicated_cdc20_not_trimmed_DM_dist(idx) , 1e2 ,'Normalization','Probability');
histogram( T.percent_underreplicated_cdc20_not_trimmed_DM_dist(~idx) , 1e2 ,'Normalization','Probability');


[f,x] = ecdf( T.percent_underreplicated_cdc20_not_trimmed_DM_dist(idx) );
plot(x,f./10,'-b','LineWidth',3)
[f,x] = ecdf( T.percent_underreplicated_cdc20_not_trimmed_DM_dist(~idx) );
plot(x,f./10,'-r','LineWidth',3)

xlim([-50 50])
xlabel('% underreplication')

[~,p] = ttest2( T.percent_underreplicated_cdc20_not_trimmed_DM_dist(idx & ~idxnan) , T.percent_underreplicated_cdc20_not_trimmed_DM_dist(~idx & ~idxnan) )


%% Best boxplot figure
fh = figure('units','centimeters','position',[5 5 4 6 ]);
idx = Q.value > prctile(Q.value,99.99) ; 
boxplot(Q.percent_underreplicated_cdc20_not_trimmed_DM_dist , idx ,'symbol','')
ylim([-9 50])
ylabel('% underreplication')
set(gca,'xtick',[])
%% Trend figure
G = (Q.value./20) ; 
G( G >= prctile(G,99.99)) = prctile(G,99.99) ;
fh = figure('units','centimeters','position',[5 5 15 6 ]);
boxplot( Q.percent_underreplicated_cdc20_not_trimmed_DM_dist , round(G)*20 ,'symbol','')
ylim([-9 50])
ylabel('% underreplication')
xlabel('G4 quadruplex score')
grid on ;

%% Top N figure
Y = Q.percent_underreplicated_cdc20_not_trimmed_DM_dist ; 
X = NaN( 2 , 1000);
for I = 1:1000
    X(1,I) = nanmedian(Y(1:I));
    X(2,I) = nanstd(Y(1:I));
end
fh = figure('units','centimeters','position',[5 5 10 6 ]);
shadedErrorBar(1:1000  ,X(1,:) , X(2,:))
set(gca,'xscale','log')
ylabel('% underreplication (median)')
xlabel('top N G4 quadruplex scoring 200nt bins')
ylim([-9 50])

%% Top N Figure v2
binsize = 12 ; 
Y = Q.percent_underreplicated_cdc20_not_trimmed_DM_dist ; 
X = NaN( 2 , 1e3);
for I = (binsize+1):(1e3+binsize)
    try
    idx = (I-binsize):(I+binsize) ;
    X(1,I-binsize) = nanmedian(Y(idx));
    X(2,I-binsize) = nanstd(Y(idx));
    catch
    end
end
fh = figure('units','centimeters','position',[5 5 10 6 ]);
shadedErrorBar( 1:size(X,2) ,  X(1,:) , X(2,:))
set(gca,'xscale','lin')
ylabel('% underreplication (median)')
xlabel(' groups of 25, G4 quadruplex scoring 200nt bins')
axis tight; 
ylim([-9 25])
xlim([0 250])

% %% playing with more analysis down here
% fn = dir('*.bed');
% for I = 1:numel(fn)
%     G4H = readtable( fn(I).name ,'Delimiter','\t','FileType','text','Format','%s%d%d%f','TreatAsEmpty' ,'.');
%     G4H = G4H( regexpcmp(G4H.Var1,'^chr[IXV]') , :);
%     G4H.Properties.VariableNames = {'chr' 'start_point' 'end_point' 'value'};
%     G4H = G4H( ~isnan(G4H.value) ,:);
%     Q = innerjoin( G4H , T , 'Key', {'chr' 'start_point' 'end_point'});
%     Q = sortrows( Q,{'chr'  'arm' 'dist_to_the_end_kb'} ,'descend');
%   
%     [c,~] = corr( Q.value ,Q.DM_diff ,'rows','complete');
%     [c1,~] = corr( Q.value(Q.dist_to_the_end_kb<100) ,Q.DM_diff(Q.dist_to_the_end_kb<100) ,'rows','complete');
%     [c2,~] = corr( Q.value(Q.dist_to_the_end_kb<50) ,Q.DM_diff(Q.dist_to_the_end_kb<50) ,'rows','complete');
% 
%     fprintf('%s\t%0.02f\t%0.02f\t%0.02f\n' , fn(I).name , c, c1, c2);
%   
% %   % I want to know if having lots of G4s along the end of an arm
% %   %  is one of the reasons for arm-to-arm differences in max %unrep
% %   %  I think this can be ansswered with cumsum() of G4 score along chr arm
% %   % but needs more thought. 
% %   Q = Q(Q.dist_to_the_end_kb<20,:);
% %     Q.ID = strcat( Q.chr , Q.arm);
% %     uchr = unique(Q.ID);
% %     Q.cumsum = NaN( height(Q),1);
% %     for uchrI = 1:numel(uchr)
% %         idx = strcmp(Q.ID,uchr{uchrI});
% %         Q.cumsum(idx) = cumsum(Q.value(idx));
% %     end
% 
% %     R = Q(Q.dist_to_the_end_kb > 0 & Q.dist_to_the_end_kb <= 15,:);
% %     G = grpstats( R , {'chr' 'arm'} , 'mean' , 'DataVars' , {'percent_unreplicated_not_trimmed_cdc20' 'value' 'cumsum'});
% %     
% %     figure; 
% %     plot( G.mean_percent_unreplicated_not_trimmed_cdc20 , G.mean_cumsum , 'ok')
% %     [c,p]=corr(G.mean_percent_unreplicated_not_trimmed_cdc20 , G.mean_cumsum ,'rows','complete')
%   % 
% 
% end
% disp('\n')
% 
