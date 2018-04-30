# Re-evaluation of G-quadruplex propensity with G4Hunter 
# Amina Bedrat  Laurent Lacroix  Jean-Louis Mergny
# Critical evidence for the biological relevance of G-quadruplexes (G4) has recently been obtained in seminal studies performed in a variety of organisms. Four-stranded G-quadruplex DNA structures are promising drug targets as these non-canonical structures appear to be involved in a number of key biological processes. Given the growing interest for G4, accurate tools to predict G-quadruplex propensity of a given DNA or RNA sequence are needed. Several algorithms such as Quadparser predict quadruplex forming propensity. However, a number of studies have established that sequences that are not detected by these tools do form G4 structures (false negatives) and that other sequences predicted to form G4 structures do not (false positives). Here we report development and testing of a radically different algorithm, G4Hunter that takes into account G-richness and G-skewness of a given sequence and gives a quadruplex propensity score as output. To validate this model, we tested it on a large dataset of 392 published sequences and experimentally evaluated quadruplex forming potential of 209 sequences using a combination of biophysical methods to assess quadruplex formation in vitro. We experimentally validated the G4Hunter algorithm on a short complete genome, that of the human mitochondria (16.6 kb), because of its relatively high GC content and GC skewness as well as the biological relevance of these quadruplexes near instability hotspots. We then applied the algorithm to genomes of a number of species, including humans, allowing us to conclude that the number of sequences capable of forming stable quadruplexes (at least in vitro) in the human genome is significantly higher, by a factor of 2–10, than previously thought.
#
#
# 
# http://dx.doi.org/10.1093/nar/gkw006
#
# https://github.com/AnimaTardeb/G4-Hunter
#

YEASTGENOME = ~/Data/Yeast/genome.fasta
dl:
	wget https://github.com/AnimaTardeb/G4-Hunter/archive/master.zip
	unzip master.zip
	mv G4-Hunter-master/* ./
	rmdir G4-Hunter-master
	rm -f master.zip

#  requires Python v2.7.5 !
w=25
rn:
	for S in 0.5 1 1.5 1.75  ; do \
		mkdir -p yeast_S$${S}_W$w ; \
		python G4Hunter_bed.py -i $(YEASTGENOME) -o yeast_S$${S}_W$w -w $w -s $${S} ;\
	done

rn2:
	for S in 0.5 1 1.5 1.75  ; do \
		bedtools map -c 5 -o sum -a ~/Projects/2017__Mendoza_ReplicationTiming/genome_windows_by_200.bed  \
			-b yeast_S$${S}_W$w/Results/genome-$w-nts | gzip > S$${S}_W$w_sum_200.bed.gz ;\
		bedtools map -c 5 -o mean -a ~/Projects/2017__Mendoza_ReplicationTiming/genome_windows_by_200.bed  \
			-b yeast_S$${S}_W$w/Results/genome-$w-nts | gzip > S$${S}_W$w_mean_200.bed.gz ;\
	done

clean:
	rm -fr yeast_S*
