% show the location of the G-quadruplexes in each chromosome, relative to the %-underreplication
% also show that more stringent thresholds results in a stronger diff
% between non-G4 & G4
DATADIR = '~/Develop/Mendoza__ReplicationEvolution/Data/';
load( [DATADIR  'DS_stat__200bp_new.mat']);
DS_FEATURES = load( [DATADIR  'DS_stat__features_new.mat']);
DS_FEATURES = DS_FEATURES.DS ; 

%% smooth out position-effects trying a few methods
DS.underrep_smoothed_1 = NaN( size(DS,1) , 1);
DS.underrep_smoothed_2 = NaN( size(DS,1) , 1);

for chrI = 1:16
    idx = DS.chr_num == chrI ;
    Y = 100*DS.percent_underreplicated_cdc20 ;
    X =  DS.middle_point_kb;
    Ys1   = smooth( X(idx) , Y(idx) ,1000,'rloess') ;
    Ys2   = smooth( X(idx) , Y(idx) ,1000,'sgolay',2) ;
    
    DS.underrep_smoothed_1(idx) = Ys1 ;
    DS.underrep_smoothed_2(idx) = Ys2 ;
    
end
DS.smooth_residual_1  = 100*DS.percent_underreplicated_cdc20 - DS.underrep_smoothed_1 ;
DS.smooth_residual_2  = 100*DS.percent_underreplicated_cdc20 - DS.underrep_smoothed_2 ;


%% plot an example of the smoothed data
chrI=4 ; 
THRESH = [99.5];
GREEN = [0.1818    0.5909    0.4000] ; 
idx = DS.chr_num == chrI ; 
fh = figure('units','centimeters','position',[5 5 15 10]); hold on ;
plot( DS.start_point(idx)./1000 , 100*DS.percent_underreplicated_cdc20(idx) ,'-','Color',GREEN)
plot( DS.start_point(idx)./1000 , DS.underrep_smoothed_1(idx) ,'-r','LineWidth',1)
plot( DS.start_point(idx)./1000 , DS.underrep_smoothed_2(idx) ,'-c','LineWidth',1)

%plot( DS.start_point(idx)./1000 , DS.smooth_residual(idx) ,'-c')

xlabel('Chromsome position (kb)')
ylabel('% under-replicated')
axis tight; 

%G quads
idx_g4_2 = DS.G4 > prctile(DS.G4,THRESH(1))  ; 
plot(DS.start_point(idx & idx_g4_2 )./1000  , 100*DS.percent_underreplicated_cdc20(idx & idx_g4_2) ,'sk','MarkerFaceColor','k')

% transposons
idx = regexpcmp(DS_FEATURES.TYPE,'transpo') & DS_FEATURES.chr_num == chrI ;
if sum(idx)>0 , plot(DS_FEATURES.middle_point_kb(idx) , -8 , 'sm') , end


title(['chr ' num2roman(chrI)]);

ylim([-9 40])

%% boxplots showing G-quads & transposon higher %underrep even after smoothing
%% G4 enrichment by random sampling
Y = DS.smooth_residual_1 ; 

THRESH = 99.9 ;
G4idx = DS.G4 > prctile(DS.G4,THRESH);
fprintf('%0.02f -> %d G-quads\n' , THRESH , sum(G4idx));
M = NaN(0);
N = 2e3 ; 
for I = 1:N
    idx = randsample( find(~G4idx) , 100 );
    M(I,1) = mean( Y(idx));
end
for I = 1:N
    idx = randsample( find(G4idx) , 100 , true );
    M(I,2) = mean( Y(idx));
end
fh = figure('units','centimeters','position',[5 5 5 8]);hold on ; 
bh = boxplot(M,'Notch','on','Symbol','') ; 
ylim([-2 prctile(M(:,2),99)])
set(gca,'xticklabel',{'random' 'G-quads'}) ; 
ylabel('Residuals (%underrep - smoothed)'); 
set(bh(5,2),'Color','k','LineWidth',2)
set(bh(5,1),'Color',GREEN,'LineWidth',2)
% Transposon enrichment
roundpos = round(DS.middle_point_kb);
TRANSidx = ismember( roundpos , round(DS_FEATURES.middle_point_kb( regexpcmp(DS_FEATURES.TYPE,'transpo'))));

fprintf('%0.02f %d transposons\n' , mean(TRANSidx) , sum(TRANSidx));
M = NaN(0);
N = 2e3 ; 
for I = 1:N
    idx = randsample( find(~TRANSidx) , 100 );
    M(I,1) = mean( Y(idx));
end
for I = 1:N
    idx = randsample( find(TRANSidx) , 100 , true );
    M(I,2) = mean( Y(idx));
end
fh = figure('units','centimeters','position',[5 5 5 8]);hold on ; 
bh = boxplot(M,'Notch','on','Symbol','') ; 
ylim([-2 prctile(M(:,2),99)])
set(gca,'xticklabel',{'random' 'tranposons'}) ; 
ylabel('Residuals (%underrep - smoothed)'); 
set(bh(5,2),'Color','m','LineWidth',2)
set(bh(5,1),'Color',GREEN,'LineWidth',2)




%% dual-y axis plot
Y = NaN(2e3 , 2);
data = 100 * DS.percent_underreplicated_cdc20 ; 
xl = linspace(95,100,nrows(Y));
for I = 1:numel(xl)
    idx =  DS.G4 > prctile(DS.G4,xl(I)) ; 
    Y(I,1) = sum(idx);
    Y(I,2) = nanmean(data(idx)) ./ nanmean(data(~idx)) ; 
 %   fprintf('%d\t%d\t%0.02f \n' , I , Y(I,1) , Y(I,2));
end

fh = figure('units','centimeters','position',[5 5 15 10 ]);
hold on; 
[AX,H1,H2] = plotyy( xl , Y(:,1) , xl , Y(:,2));
set(H1,'LineWidth',3)
set(H2,'LineWidth',3)

grid on ;
set(gca,'ytick',[1 10 20 50 100 200 500 1e3 2e3 5e3])
set(AX(1),'yscale','log')
ylabel(AX(1)','# of G-quadruplexes');
ylabel(AX(2)','%underrep (G-quad / non-G-quad)');
xlabel('threshold value (% of G4-Hunter scores in the genome)')
set(AX(2),'ytick',1:50)


%% dual-y axis plot for SNPs shows that it's just a position thing
Y = NaN(2e3 , 2);
data = 100 * DS.percent_underreplicated_cdc20 ; 
data = DS.smooth_residual_1 ;
xl = linspace(50,100,nrows(Y));
for I = 1:numel(xl)
    idx =  DS.freq_SNP > prctile(DS.freq_SNP,xl(I)) ; 
    Y(I,1) = sum(idx);
    Y(I,2) = nanmean(data(idx)) ; 
    Y(I,3) = nanmean(data(~idx)) ; 
end

fh = figure('units','centimeters','position',[5 5 15 10 ]);
hold on; 
[AX,H1,H2] = plotyy( xl , Y(:,1) , xl , Y(:,2));
set(H1,'LineWidth',3)
set(H2,'LineWidth',3)
grid on ;
set(AX(1),'ytick',[1 10 20 50 100 200 500 1e3 2e3 5e3 1e4 2e4])
set(AX(1),'yscale','log')
ylabel(AX(1)','# of loci');
ylabel(AX(2)','%underrep - smoothed');
xlabel('threshold value (% high SNP loci)')
%set(AX(2),'ytick',1:50)
line( AX(2) , xlim , [0 0],'LineStyle','--','Color',[.7 .7 .7])
set(AX(2),'ylim',[-0.2 2.2])
set(AX(2),'ytick',0:5)


% 
% %% make new table w/DS & DS_FEATURES
% DS.middle_point_kb = round(DS.middle_point_kb);
% DS_FEATURES.middle_point_kb = round(DS_FEATURES.middle_point_kb);
% 
% vn = {'smooth_residual' 'Trep_spline' 'dist_to_ARS' 'GC' 'G4' 'IS' 'underrep_smoothed' 'subtelomeric_bool' 'max_PROseq' ...
%     'median_freq_N_chr_wo_gaps' 'median_freq_N_chr_w_gaps' 'freq_SNP' 'freq_indel' } ; 
% G = grpstats( dataset2table(DS) , {'chr_num' 'middle_point_kb'}, {'nanmean' 'max'} , 'DataVars'  , vn);
% G.transposon = false(height(G),1);
% for I = 1:height(G)
%     G.transposon(I) = any(regexpcmp( DS_FEATURES.TYPE( DS_FEATURES.chr_num == G.chr_num(I) & DS_FEATURES.middle_point_kb == G.middle_point_kb(I)) , 'transp')) ; 
% end
% %%
% vn = {'smooth_residual' 'Trep_spline' 'dist_to_ARS' 'GC' 'G4' 'max_PROseq' } ; 
% vn = [ strcat('nanmean_' , vn(2:end)) 'transposon'] ;
% Y = G.nanmean_smooth_residual  ; 
% X = table2array(G( : , vn ))  ; 
% [P,T,STATS,TERMS]=anovan( Y , X ,'Continuous',1:numel(vn)-1 , 'VarNames' , vn );
% 
