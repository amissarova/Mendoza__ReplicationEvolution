%% look at coverage for all cerevisae gDNA sequencing data
% needs depth.tgz from the cluster
% uncompress into a folder and run this script
% the script does some filtering, makes some figures, and saves a .mat file w/the filtered experiments. 


%% load coverage data for all experiments
% download & uncompress depth.tgz as generated from the 'depth' target in 
%    /homes/users/lcarey/single_cell_behavior/Projects/2017__Mendoza_ReplicationTiming/AllYeastGDNAseqs
% eg: 
%   scp marvin.s.upf.edu:Projects/2017__Mendoza_ReplicationTiming/AllYeastGDNAseqs/depth.tgz ./ ; tar -czvf depth.tgz

DATADIR = '~/Desktop/AllYeastGDNAseqs/' ; 
EXPS = readtable([ DATADIR 'exps.tab' ] ,'FileType','text','ReadVariableNames',false);
DEPTH = dlmread( [ DATADIR 'depth.tab' ] )';
LOCS = readtable([ DATADIR 'regions.tab' ],'FileType','text','ReadVariableNames',false);

LOCS.Properties.VariableNames = {'chr' 'start' 'stop'};

%remove experiments w/zero coverage at all positions (incomplete files)
%EXPS  =  EXPS( ~all(DEPTH==0,2) ,:);
%DEPTH = DEPTH( ~all(DEPTH==0,2) ,:);

% remove non-gDNA seqs
DEPTH = DEPTH( regexpcmp( LOCS.chr , '^chr[IVX]') ,:); 
LOCS = LOCS( regexpcmp( LOCS.chr , '^chr[IVX]') ,:); 

%% remove experiments w/uneven coverage
%   arbitrary filters, but the results seem rather insensitive to reasonable values

% filter to get experiments w/reasonable coverage
THRESH_MEAN_VS_MEDIAN = 1.05 ; % 1.05 very stringent , 1.5 loose
%THRESH_MEAN_VS_MEDIAN = 10.5 ; % 1.05 very stringent , 1.5 loose

THRESH_MAX_BINS_ZERO = numel(LOCS.chr)/10 ; % discard experiment if 10% is zero


normalized_depth = DEPTH ./ repmat( sum(DEPTH) , size(DEPTH,1) , 1);
medians = median(normalized_depth,1) ; 
means   = mean(normalized_depth,1) ; 

idx_reasonable_distribution = (medians > (means./THRESH_MEAN_VS_MEDIAN)) & (medians < (means.*THRESH_MEAN_VS_MEDIAN)) ;  % stringent threshold determined by visual inspection of imagesc()

idx_not_missing_too_much = sum(normalized_depth==0) < THRESH_MAX_BINS_ZERO ; % no one that passes the 'reasonable distribution' threshold is missing > 200kb
idx = idx_not_missing_too_much & idx_reasonable_distribution ; 

EXPS_passed_thresh = EXPS.Var1(idx, :) ; 
normalized_depth = normalized_depth( : , idx ) ;
figure; imagesc( normalized_depth , [0 prctile(normalized_depth(:),95)] )

%% find the % of samples for which each genomic window is > 1,2,3 fold above/below the median coverage for that position
uchr = unique(LOCS.chr);
R = table();
for I = 1:numel(uchr)
    Q = table();
    this_chr_idx = strcmp(LOCS.chr,uchr{I});
    Q.chrlocs = LOCS.start(this_chr_idx) ; 
    Q.chr = repmat(uchr(I) , height(Q) , 1);
    Q.kb_from_end = min( [ (Q.chrlocs ./ 1000) ( max(Q.chrlocs) - Q.chrlocs ) ./ 1000 ] , [] , 2);

    X = normalized_depth(this_chr_idx,:);
    median_CN_each_exp = median(X,1);

    X_norm_div_per_exp_median = X ./ repmat(median_CN_each_exp,size(X,1),1) ; 

    median_CN_per_pos = median( X_norm_div_per_exp_median , 2 );

    CN_l2 = log2( X_norm_div_per_exp_median ./ repmat( median_CN_per_pos , 1 , size(X,2)) ) ;

    Q.pct_above_1 = 100 * mean(CN_l2>1  , 2)  ;
    Q.pct_above_2 = 100 * mean(CN_l2>2  , 2)  ;
    Q.pct_above_3 = 100 * mean(CN_l2>3  , 2)  ;
   
    
    Q.pct_below_1 = 100 * mean(CN_l2<-1  , 2)  ;
    Q.pct_below_2 = 100 * mean(CN_l2<-2  , 2)  ;
    Q.pct_below_3 = 100 * mean(CN_l2<-3  , 2)  ;
   

    R = vertcat(Q,R);
end
R = sortrows(R , 'kb_from_end','ascend');
%% figure
% not sure if mean or median makes the most sense

xl = unique(R.kb_from_end( R.kb_from_end < 100));
data = NaN( 1e3 , numel(xl));
for I = 1:size(data,2)
    data( : , I) = bootstrp( size(data,1) , @median , R.pct_above_1(R.kb_from_end==xl(I)));
end


fh = figure('units','centimeters','position',[5 5 8 8 ]);
h = shadedErrorBar( xl , median(data) , std(data)) ; 
set(h.mainLine,'LineWidth',2) ;
xlabel('KB from the telomere')
ylabel('% of samples w/CN > 2')
ylim([0 max(ylim)])
ylim([0 10])
set(gca,'ytick',0:10)
grid on ;
%%


% %%
% THRESH = [ 0.15 1 5] ; %Initial CN thresh ; mean CN thresh ; run-lenth thresh
% figure; hold on ;
% for I = 1
%     l2_CN_this_exp = CN_l2(:,I) ; 
%     [start_idx,end_idx,run_len] = consecutive_ones( l2_CN' > THRESH(1) ) ; 
%     cnruns = arrayfun( @(I)median(l2_CN(start_idx(I):end_idx(I))) , 1:numel(start_idx)) ;
%     idx = cnruns > THRESH(2) & run_len >= THRESH(3) ;
%     start_idx = start_idx(idx);
%     end_idx = end_idx(idx);
%     run_len = run_len(idx);
% 
%     plot( chrlocs , l2_CN , '.');
%     plot( chrlocs(start_idx:end_idx) , l2_CN(start_idx:end_idx) , '-o');
%     
% end
% ylim([-1 2])
% %%
% R = table();
% warning('off','MATLAB:table:RowsAddedExistingVars') ; 
% THRESH = prctile(normalized_depth(:),1) ; 
% chrs = unique(LOCS.chr);
% c=0;
% for chri = 1:numel(chrs)
%     idx = find(strcmp(LOCS.chr,chrs{chri}));
%     for I = 1:1000
%         missing_this_much = find(all(normalized_depth( idx(1:I) , : )  <= THRESH, 1)) ; 
%         if any(missing_this_much)
%             c=c+1;
%             R.chr{c} = chrs{chri} ;
%             R.topos(c) = LOCS.stop(idx(I)) ;
%             R.nmissing(c) =  numel(missing_this_much);
%             R.whomissing{c} =  EXPS_passed_thresh(missing_this_much) ;
%             fprintf( '%s\t0-%d\t%d\t' , chrs{chri} , LOCS.stop(idx(I)) ,  numel(missing_this_much) );
%             for exp = missing_this_much
%                 fprintf('%s\t' , EXPS_passed_thresh{exp})
%             end
%             fprintf('\n');
%             
%         end
%     end
% end
% 
% % plot results
% figure;
% hold on; 
% clrs = parula(numel(chrs));
% BINSIZE = 200 ; 
% markers = 'o+svo+svo+svo+sv';
% for I = 1:numel(chrs)
%     idx = find(strcmp(R.chr,chrs{I}));
%     if ~isempty(idx)
%     X = [ R.topos(idx)' max(R.topos(idx))+BINSIZE]; 
%     Y =  [ R.nmissing(idx)' 0] ; 
%     Y = (Y ./ ncols(normalized_depth) ) * 100  ; % % of strains
%     plot( X ./ 1000 , Y,'-' ,'Marker',markers(I), 'Color',clrs(I,:),'DisplayName',chrs{I} ,'LineWidth',3);
%     end
% end
% ylabel('% of strains missing this')
% xlabel('Kb from left end')
% set(gca,'xscale','log')
% axis tight;
% set(gca,'xtick',[0.1 0.2 0.5 1 2 5 10])
% legend('location','ne')
% 
% %% save results
% save('AllYeastGDNAseqs.mat' , 'R' , 'normalized_depth'  , 'EXPS_passed_thresh' , 'LOCS') ;
