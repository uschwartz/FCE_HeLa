## input

+ takes a tab-delimited table containing the GC normalization factors derived from loess regression 

|chrom  |start  |end  |normFactor |
|---	|---	|---	|---	|
|chr20   	|1	   	|250   	|3   	|
|chr20   	|251	   	|500  	|4.5   	|
|chr20   	|751	  	|750   	|1.2   	|

...

+ takes a wig file (e. g. DANPOS2 output)

|fixedStep chrom=chr20 start=1 step=10 span=10 |
|---	|
|3 |
|3 |
|3 |
|6 |
|6 |

...


## output

normalized wig file

|fixedStep chrom=chr20 start=1 step=10 span=10 |
|---	|
|1 |
|1 |
|1 |
|2 |
|2 |

...
