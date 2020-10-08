## input

+ takes a tab-delimited table containing the GC normalization factors derived from loess regression 

chrom start	end	normFactor
chr20	1	250	0
chr20	251	500	0
chr20	501	750	0
chr20	751	1000	0
...

+ takes a wig file (e. g. DANPOS2 output)

fixedStep chrom=chr20 start=1 step=10 span=10
0
0
0
0
0
0
0
...


## output

normalized wig file

fixedStep chrom=chr20 start=1 step=10 span=10
0
0
0
0
0
0
0
...
