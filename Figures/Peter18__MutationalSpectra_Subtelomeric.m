%% load data
DATADIR = '~/Data/Peter18/' ;
T = readtable([ DATADIR 'nocomma.tab' ] ,'FileType','text','Delimiter','\t','Format','%d%d%s%s');
T.Properties.VariableNames = {'chr_num' 'pos' 'ref' 'alt'};
%% Define Subtelomeric regions
SUBTEL = 25*1000 ; 
T.Subtel = zeros(height(T),1);
T.Subtel( T.pos < SUBTEL) = 1 ; 
T.Subtel( T.pos > (max(T.pos)-SUBTEL)) = 1 ; 

%% keep only SNPs
T = T( cellfun(@length,T.alt)==1 , :);
T = T( cellfun(@length,T.ref)==1 , :);


%% summary stats
G=grpstats(T,{'chr_num' 'ref' 'alt' 'Subtel'}, 'sum','DataVars','chr_num');
G.sum_chr_num = []; 


%% compute stats for each mutation & chr
G.FE_p = NaN(height(G),1);
G.FE_or = NaN(height(G),1);
for I = 1:2:height(G)
    st = sum(G.GroupCount( G.Subtel==1 & G.chr_num==G.chr_num(I)));
    nonst = sum(G.GroupCount( G.Subtel==0 & G.chr_num==G.chr_num(I)));
    [~,p,or] = fishertest( [ nonst st ; G.GroupCount(I) G.GroupCount(I+1)] );
   G.FE_p(I) = log10(p);
   G.FE_or(I) = or.OddsRatio ;
end


G2 = grpstats( G( ~isnan(G.FE_p) ,:) , {'ref' 'alt'},{'mean' 'std'},'DataVars',{'FE_p' 'FE_or'})

% %     %%
% % all_subel__nonsubtel = [sum(G.GroupCount(G.Subtel25==1)) sum(G.GroupCount(G.Subtel25==0))];
% % figure;
% % hold on ;
% % for I = 1:2:height(G)
% %     non_subtel__count =  G.GroupCount(I) ;
% %     subtel__count =  G.GroupCount(I+1) ;
% %     txt = sprintf( '%s->%s' , G.Var3{I}  , G.Alt1{I} );
% %     [~,p,or] = fishertest( [ all_subel__nonsubtel ; subtel__count non_subtel__count]) ;
% %     plot( non_subtel__count , subtel__count , '.k');
% %     if p < 0.05
% %         text( non_subtel__count , subtel__count , txt);
% %     end
% % end
% % xlabel('non-subtel count ')
% % ylabel('subtel count')
% % %%
% % lr = log2(all_subel__nonsubtel(1)/all_subel__nonsubtel(2)); 
% % for I = 1:2:height(G)
% %     non_subtel__count =  G.GroupCount(I) ;
% %     subtel__count =  G.GroupCount(I+1) ;
% %     [~,p,or] = fishertest( [ all_subel__nonsubtel ; subtel__count non_subtel__count]) ;
% %     txt = sprintf( '%s->%s' , G.Var3{I}  , G.Alt1{I} );
% %     fprintf( '%s\t%0.02f\t%0.03f\n' , txt, log2(subtel__count./non_subtel__count) - lr  , p );
% % end
% % 
% % %%
