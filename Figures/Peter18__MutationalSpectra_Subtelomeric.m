T = readtable('x2.tab','FileType','text');
T.Subtel50 = zeros(height(T),1);
T.Subtel25 = zeros(height(T),1);
T.Subtel50( T.Var2 < 50000) = 1 ; 
T.Subtel25( T.Var2 < 25000) = 1 ; 
T.Subtel50( T.Var2 < (max(T.Var2)-50000)) = 1 ; 
T.Subtel25( T.Var2 < (max(T.Var2)-25000)) = 1 ; 
T.Alt1 = regexprep(T.Var4,',.*','');
%%
AllMuts=grpstats(T,{'Subtel25'},'sum','DataVars','Var1');
G=grpstats(T,{'Var3' 'Alt1' 'Subtel25'},'sum','DataVars','Var1');
G = G( cellfun(@length,G.Alt1)==1 , :);
G = G( cellfun(@length,G.Var3)==1 , :);
%%
all_subel__nonsubtel = [sum(G.GroupCount(G.Subtel25==1)) sum(G.GroupCount(G.Subtel25==0))];
figure;
hold on ;
for I = 1:2:height(G)
    non_subtel__count =  G.GroupCount(I) ;
    subtel__count =  G.GroupCount(I+1) ;
    txt = sprintf( '%s->%s' , G.Var3{I}  , G.Alt1{I} );
    [~,p,or] = fishertest( [ all_subel__nonsubtel ; subtel__count non_subtel__count]) ;
    plot( non_subtel__count , subtel__count , '.k');
    if p < 0.05
        text( non_subtel__count , subtel__count , txt);
    end
end
xlabel('non-subtel count ')
ylabel('subtel count')
%%
lr = log2(all_subel__nonsubtel(1)/all_subel__nonsubtel(2)); 
for I = 1:2:height(G)
    non_subtel__count =  G.GroupCount(I) ;
    subtel__count =  G.GroupCount(I+1) ;
    [~,p,or] = fishertest( [ all_subel__nonsubtel ; subtel__count non_subtel__count]) ;
    txt = sprintf( '%s->%s' , G.Var3{I}  , G.Alt1{I} );
    fprintf( '%s\t%0.02f\t%0.03f\n' , txt, log2(subtel__count./non_subtel__count) - lr  , p );
end

%%
T = readtable('~/Downloads/nocomma.tab','FileType','text);