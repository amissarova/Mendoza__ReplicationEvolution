tobed:
	modify_column.pl -m 1000 -c 1 < Alvino07.txt | cut.pl -f 1,2,2,3 |grep -v PctRep_40 > Alvino07.bed

200:
	bedtools  map -a genome_200bpwindows_GC.bed -b Raghuraman01.bed -c 4 -o mean -null SKIPME | grep -v SKIPME|cut -f 1,2,5 > Raghuraman01_200.txt
	bedtools  map -a genome_200bpwindows_GC.bed -b Koren10_WT.bed -c 4 -o mean -null SKIPME | grep -v SKIPME|cut -f 1,2,5 > Koren10_WT_200.txt
	bedtools  map -a genome_200bpwindows_GC.bed -b Alvino07.bed -c 4 -o mean -null SKIPME | grep -v SKIPME|cut -f 1,2,5 > Alvino07_200.txt
