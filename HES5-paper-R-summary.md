HES5 methylation in prostate cancer: **code** and **plots**
========================================================

This document contains descriptions and code required to reproduce the figures from our paper. The data sets required to reproduce these figures are bundled together with the R Markdown document which was used to create this document.  The R Markdown document can be parsed using 'knitr' package in R to reproduce this document and all the figures from our paper, if you so wish.

___________
___________
**Figure 1**
-------------------------

**Code and plot for Figure 1b:**



```r
load("ComplexMeth.RData")
ls()
```

```
##  [1] "B006"       "B007"       "B008"       "B054"       "BBlood"    
##  [6] "GMAnno"     "hm450"      "M006"       "M007"       "M008"      
## [11] "M054"       "MBlood"     "PriceAnno"  "PriceAnno2" "targets"
```

```r
library(GenomicRanges)
```

```
## Loading required package: BiocGenerics
## Loading required package: parallel
## 
## Attaching package: 'BiocGenerics'
## 
## The following objects are masked from 'package:parallel':
## 
##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
##     clusterExport, clusterMap, parApply, parCapply, parLapply,
##     parLapplyLB, parRapply, parSapply, parSapplyLB
## 
## The following object is masked from 'package:stats':
## 
##     xtabs
## 
## The following objects are masked from 'package:base':
## 
##     anyDuplicated, append, as.data.frame, as.vector, cbind,
##     colnames, do.call, duplicated, eval, evalq, Filter, Find, get,
##     intersect, is.unsorted, lapply, Map, mapply, match, mget,
##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
##     rbind, Reduce, rep.int, rownames, sapply, setdiff, sort,
##     table, tapply, union, unique, unlist
## 
## Loading required package: IRanges
## Loading required package: GenomeInfoDb
```

```r
gene.list <- rbind(c("HES5","chr1", 2461684, 2460184,"NM_001010926"), 
                   c("APC", "chr5", 112073556, 112181936, "NM_001127510"), 
                   c("RARB", "chr3",25469834, 25639422 , "NM_000965" ), 
                   c("TACC2", "chr10",123922941, 124014060, "NM_001291878"), 
                   c("ITGB2", "chr21", 46340965, 46305868, "NM_000211" ),
                   c("DGKZ", "chr11", 46383145, 46402104, "NM_001105540"), 
                   c("C5orf4", "chr5" , 154230213, 154198052, "NM_032385"), 
                   c("mir10b", "chr2",177015031, 177015140, "NR_029609"))

for(i in 1:dim(gene.list)[1]){
  
  GENE <- gene.list[i,1]
  CHR <- gene.list[i,2]
  START <- as.numeric(gene.list[i,3])-200
  END <- as.numeric(gene.list[i,3])+200
  
  assign(paste(GENE, "B006.TNmed", sep="."), median(B006[which(PriceAnno$MAPINFO.1>=START & PriceAnno$MAPINFO.1<=END & seqnames(hm450)==CHR),1:6]/B006[which(PriceAnno$MAPINFO.1>=START & PriceAnno$MAPINFO.1<=END & seqnames(hm450)==CHR),7]))

  assign(paste(GENE, "B007.TNmed", sep="."), median(B007[which(PriceAnno$MAPINFO.1>=START & PriceAnno$MAPINFO.1<=END & seqnames(hm450)==CHR),1:7]/B007[which(PriceAnno$MAPINFO.1>=START & PriceAnno$MAPINFO.1<=END & seqnames(hm450)==CHR),8]))
  
  assign(paste(GENE, "B008.TNmed", sep="."), median(B008[which(PriceAnno$MAPINFO.1>=START & PriceAnno$MAPINFO.1<=END & seqnames(hm450)==CHR),1:3]/B008[which(PriceAnno$MAPINFO.1>=START & PriceAnno$MAPINFO.1<=END & seqnames(hm450)==CHR),4]))
  
  assign(paste(GENE, "B054.TNmed", sep="."), median(B054[which(PriceAnno$MAPINFO.1>=START & PriceAnno$MAPINFO.1<=END & seqnames(hm450)==CHR),1]/B054[which(PriceAnno$MAPINFO.1>=START & PriceAnno$MAPINFO.1<=END & seqnames(hm450)==CHR),2]))
}

B006.TN.medians.for.heatmap <- rbind(HES5.B006.TNmed, APC.B006.TNmed, RARB.B006.TNmed, TACC2.B006.TNmed, ITGB2.B006.TNmed, DGKZ.B006.TNmed, C5orf4.B006.TNmed, mir10b.B006.TNmed)
B007.TN.medians.for.heatmap <- rbind(HES5.B007.TNmed, APC.B007.TNmed, RARB.B007.TNmed, TACC2.B007.TNmed, ITGB2.B007.TNmed, DGKZ.B007.TNmed, C5orf4.B007.TNmed, mir10b.B007.TNmed)
B008.TN.medians.for.heatmap <- rbind(HES5.B008.TNmed, APC.B008.TNmed, RARB.B008.TNmed, TACC2.B008.TNmed, ITGB2.B008.TNmed, DGKZ.B008.TNmed, C5orf4.B008.TNmed, mir10b.B008.TNmed)
B054.TN.medians.for.heatmap <- rbind(HES5.B054.TNmed, APC.B054.TNmed, RARB.B054.TNmed, TACC2.B054.TNmed, ITGB2.B054.TNmed, DGKZ.B054.TNmed, C5orf4.B054.TNmed, mir10b.B054.TNmed)

B006.7.8.54.TN.medians.for.heatmap <- cbind(B006.TN.medians.for.heatmap, B007.TN.medians.for.heatmap, B008.TN.medians.for.heatmap, B054.TN.medians.for.heatmap)
colnames(B006.7.8.54.TN.medians.for.heatmap) <- c("B006", "B007", "B008", "B054")

heatmap(B006.7.8.54.TN.medians.for.heatmap, scale='col', Colv=NA, margins=c(5,10), labRow=c(gene.list[,1]))
```

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1.png) 



**Code and plot for Figure 1c:**


```r
col6ii <- c(rep("red",6), "black", "dark grey")
col7ii <- c(rep("red",7), "black", "dark grey")
col8ii <- c(rep("red",3), "black", "dark grey")
col54ii <- c("red", "black")

GENE <- "HES5"
CHR <- "chr"
START <- 2461684-200
END <- START+400

par(mfrow=c(1,4), mar=c(6,3.75,2,0), mgp=c(2.5,1,0))
boxplot(cbind(B006[which(PriceAnno$MAPINFO.1>=START & PriceAnno$MAPINFO.1<=END),],BBlood[which(PriceAnno$MAPINFO.1>=START & PriceAnno$MAPINFO.1<=END),1]), las=2, main="Case-006", names=c(colnames(B006[which(PriceAnno$MAPINFO.1>=START & PriceAnno$MAPINFO.1<=END),]), "Blood"), ylab=paste(GENE,"promoter methylation (B-value)"), border=col6ii, ylim=c(0,1))

par(mar=c(6,3.5,2,0))
boxplot(cbind(B007[which(PriceAnno$MAPINFO.1>=START & PriceAnno$MAPINFO.1<=END),],BBlood[which(PriceAnno$MAPINFO.1>=START & PriceAnno$MAPINFO.1<=END),2]), las=2, main="Case-007", names=c(colnames(B007[which(PriceAnno$MAPINFO.1>=START & PriceAnno$MAPINFO.1<=END),]), "Blood"), border=col7ii, ylim=c(0,1))

par(mar=c(6,3.5,2,2))
boxplot(cbind(B008[which(PriceAnno$MAPINFO.1>=START & PriceAnno$MAPINFO.1<=END),],BBlood[which(PriceAnno$MAPINFO.1>=START & PriceAnno$MAPINFO.1<=END),3]), las=2, main="Case-008", names=c(colnames(B008[which(PriceAnno$MAPINFO.1>=START & PriceAnno$MAPINFO.1<=END),]), "Blood"), border=col8ii, ylim=c(0,1))

par(mar=c(6,1.5,2,7.5))
boxplot(B054[which(PriceAnno$MAPINFO.1>=START & PriceAnno$MAPINFO.1<=END),], las=2, names=c("T1", "Benign"), main="Case-054", border=col54ii, ylim=c(0,1))
```

![plot of chunk unnamed-chunk-2](figure/unnamed-chunk-2.png) 

___________

**Code and plots for Figure 1d-e:**


```r
library(parallel)
library(BiocGenerics)
library(IRanges)
library(XVector)
library(Biostrings)
library(XVector)
library(GenomicRanges)

col6 <- c("red", "red", "red", "dark blue", "light green","orange", "black", "dark grey")
col7 <- c("red", "dark blue", "dark blue", "dark blue", "light green", "orange", "pink", "black", "dark grey")
col8 <- c("red", "dark blue", "light green", "black", "dark grey")
col54 <- c("red", "black")

Fig.1c.1d <- rbind(c("HES5","chr1" ,2460184, 2461684 ,"NM_001010926", "-"),
c("GSTP1","chr11" ,67351066, 67354124 ,"NM_000852", "+"))

for(i in 1:2){
GENE.bucket <- Fig.1c.1d[i,]
WIN <- 5000
GENE <- GENE.bucket[1]
CHR <- GENE.bucket[2]
START <- as.numeric(GENE.bucket[3])
END <- as.numeric(GENE.bucket[4])
NAME <- GENE.bucket[5]
GENE.PROBE.ROWS <- which(seqnames(hm450)==CHR & start(hm450)>START-WIN & end(hm450)<START+WIN)

par(mar=c(7,6,4,2), mgp=c(5,0.7,0), mfrow=c(1,1))
plot(x=(PriceAnno[GENE.PROBE.ROWS,4]), y=log(B006[GENE.PROBE.ROWS,1]/B006[GENE.PROBE.ROWS,7]), type="b", col="red", main=paste(GENE, " locus"), cex=0.3, xlab=paste("Genomic coordinates of ",GENE," probes"), ylab="Methylation change\nlog(Tumour/benign)", las=2, ylim=c(-2,4))

points(x=(PriceAnno[GENE.PROBE.ROWS,4]), y=log(B006[GENE.PROBE.ROWS,2]/B006[GENE.PROBE.ROWS,7]), type="b", col="red", cex=0.3)

points(x=(PriceAnno[GENE.PROBE.ROWS,4]), y=log(B006[GENE.PROBE.ROWS,3]/B006[GENE.PROBE.ROWS,7]), type="b", col="red", cex=0.3)

points(x=(PriceAnno[GENE.PROBE.ROWS,4]), y=log(B006[GENE.PROBE.ROWS,4]/B006[GENE.PROBE.ROWS,7]), type="b", col="light blue", cex=0.3)

points(x=(PriceAnno[GENE.PROBE.ROWS,4]), y=log(B006[GENE.PROBE.ROWS,5]/B006[GENE.PROBE.ROWS,7]), type="b", col="light green", cex=0.3)

points(x=(PriceAnno[GENE.PROBE.ROWS,4]), y=log(B006[GENE.PROBE.ROWS,6]/B006[GENE.PROBE.ROWS,7]), type="b", col="orange", cex=0.3)

points(x=(PriceAnno[GENE.PROBE.ROWS,4]), y=log(B006[GENE.PROBE.ROWS,7]/B006[GENE.PROBE.ROWS,7]), type="b", col="black", cex=0.3)

lines(x=c(START,END), y=c(-1.5, -1.5))
text(x=END-1000, y=-1.9, labels=NAME)
legend(x=START-3500, y=4, legend=c(colnames(B006)), lty=1, col=col6, cex=0.6)

}
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-31.png) ![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-32.png) 
___________
___________

**Figure 2**
-------------------------

**Code and plots for Figure 2a:**

Fastq files from targeted HES5 promoter bisulfide sequencing of 39 matched Tumour-Normal pairs were demultiplexed and reads were mapped to a custom HES5 reference using Bismark:

<font color=#A4A4A4>
>#from folder containing .fasta reference sequence (or genome folder) [run once to build BS-reference]
>
>bismark_v0.12.2/bismark_genome_preparation .
>
>#from folder containing BS fastq files (outputs .sam and report.txt):
>
>../bismark_v0.12.2/bismark -n 1 -l 50 ../ -1 SLX-8224.N701_N507.000000000-A8R33.s_1.r_1.fq -2 SLX-8224.N701_N507.000000000-A8R33.s_1.r_2.fq 
>
>#from same folder (outputs CpG, CHH, CHG reports):
>
>../bismark_v0.12.2/bismark_methylation_extractor -s --comprehensive SLX-8224.N701_N507.000000000-A8R33.s_1.r_1.fq_bismark_pe.sam 
>
</font>


To process all 98 demultiplexed samples (39xTumour + 39xBenign) three simple perl scripts were written (*bismarkwrapper.pl*, *methextractwrapper.pl*, *percentMeth.pl*):

***bismarkwrapper.pl***
<p> <i>

#!/usr/bin/perl

use strict;
use warnings;

my $lof = "fqList\.txt";
my @files = ();
my $i;

open( INFILE, "<$lof")
        or die( "File $lof not found\.\n");
@files = <INFILE>;
close( INFILE );

 
for( $i = 0; $i <=  $#files; $i +=2 ){ 
my $file1="";
my $file2="";
  	$file1 = $files[$i];
		$file2 = $files[$i+1];
		chomp($file1);
		chomp($file2);
		print( "Processing $file1 and $file2\n" );
		system("\.\./bismark_v0\.12\.2/bismark -n 1 -l 50 \.\./ -1 $file1 -2 $file2");
		
}

</i></p>

Run like this:

<font color=#A4A4A4>
>ls fastqData/*.gz > fqList.txt
>
>perl ./bismarkwrapper.pl
>
>#outputs all .sam and report.txt files
</font>

***methextractwrapper.pl***

<p><i>
#!/usr/bin/perl
use strict;
use warnings;

my $lof = "samList\.txt";
my @files = ();
my $i = "";

open( INFILE, "<$lof")
  or die( "Couldn't opne file $lof" );
@files = <INFILE>;
close( INFILE );

for $i (0 .. $#files){
	chomp($files[$i]);
	print( "Extracting methylation calls for $files[$i]\.\n" );
	system("\.\./bismark_v0\.12\.2/bismark_methylation_extractor -s  --comprehensive $files[$i]");
	}
</i></p>

Run like this:

<font color=#A4A4A4>
>ls *report.txt > reportList.txt
>
>ls *.sam > samList.txt
>
>perl ./methextractwrapper.pl
>
>#outputs CpG, CHG and CHH reports
</font>



***percentMeth.pl***

<p><i>
#!/usr/bin/perl
use warnings;
use strict;

my $file = "percent-methylation\.txt" ;
my $file2 = "formatted-percent-methylation\.txt";
my @data = ();
my @line = ();
my $i ;

open( INFILE, "<", $file )
  or die( "Couldn't open file $file\.\n" ) ;
@data = <INFILE>;
close( INFILE );

open( OUTFILE, ">", $file2 )
	or die( "Couldn't open file $file2\.\n");

for $i( 0 .. $#data ){
	my @line = split(/\t/, $data[$i]);
	my @names = split(/\./, $line[0]);
	chomp($names[1]);
	chomp($line[1]);
	chop($line[1]);
	print( OUTFILE "$names[1]\t$line[1]\n");
	}

close( OUTFILE );
</i></p>

Run like this:

<font color=#A4A4A4>
>#get %meth data out:
>
>grep 'C methylated in CpG context:' * > percent-methylation.txt
>
>#after grep-ing lines out for %meth fix format for R plots and decode:
>
>perl ./percentMeth.pl
</font>


Extract unique coverage:

<font color=#A4A4A4>
>grep 'Number of paired-end alignments with a unique best hit:' * > percent-methylation.txt
</font>


Finally, select mCpG calls from Bismark CpG files using awk:

<font color=#A4A4A4>
>awk '{OFS="\t"; print $4, $5}' CpG_context_SLX-8224.N702_N507.000000000-A8R33.s_1.r_1.fq.gz_bismark_pe.txt > N702_N507_CpG_calls.txt
>
>awk '{OFS="\t"; print $4, $5}' CpG_context_SLX-8224.N707_N506.000000000-A8R33.s_1.r_1.fq.gz_bismark_pe.txt > N707_N506_CpG_calls.txt
</font>

Then read into R and visualize using mosaic plot:


```r
load("hes5Validation2.Rdata")


#read in:
#TB09.1008_T
#N702_N507.CpG <- read.table("N702_N507_CpG_calls.txt", sep="\t", skip=1)
#TB09.1008_B
#N707_N506.CpG <- read.table("N707_N506_CpG_calls.txt", sep="\t", skip=1)

#plot CpG methylation contingency table (mosaic-plot:)
plot(table(N707_N506.CpG), col=c("blue", "red"), border=c("blue", "red"),  main="TB09.1008 Benign\nHES5 promotor bisulphite sequencing", xlab="CpG positions", ylab="Proportion methylated / unmethylated")
```

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-41.png) 

```r
plot(table(N702_N507.CpG), col=c("blue", "red"), border=c("blue", "red"),  main="TB09.1008 Tumour\nHES5 promotor bisulphite sequencing", xlab="CpG positions", ylab="Proportion methylated / unmethylated")
```

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-42.png) 
___________

**Code and plots for Figure 2b:**


```r
#plot-xy:
plot(x=full.data.b.core[order(full.data.b.core[,3]),2]/100, y=full.data.t.core[order(full.data.t.core[,3]),2]/100, xlim=c(0,1), ylim=c(0,1), cex=0.5, pch=19, xlab="Benign methylation", ylab="Tumour methylation", main="HES5 promoter methylation\npaired T-N bisulphite sequencing")

abline(0,1, lty=2)
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-5.png) 

___________

**Code and plots for Figure 2c:**


```r
samples
```

```
##  [1] "TB09.0817" "TB09.1008" "TB09.1083" "TB09.1254" "TB09.1365"
##  [6] "TB09.1402" "TB09.1424" "TB10.1062" "TB10.1284" "TB10.1416"
## [11] "TB11.0340" "TB11.0387" "TB11.0598" "TB11.0675" "TB11.0695"
## [16] "TB11.1035" "TB11.1083" "TB11.1437" "TB11.1447" "TB11.1693"
## [21] "TB11.1766" "TB11.1889" "TB11.1912" "TB11.2182" "TB11.2279"
## [26] "TB11.2442" "TB12.0135" "TB12.0144" "TB12.0170" "TB12.0322"
## [31] "TB12.0392" "TB12.0518" "TB12.1231" "TB12.1447" "TB12.1505"
## [36] "TB12.1843" "TB12.1852" "TB12.2315" "TB12.2336"
```

```r
hes5.wilcox <- data.frame(samples)
hes5.wilcox[,2] <- c(rep(1, 39))
names(hes5.wilcox) <- c("samples", "p-values")

for(i in 1:length(samples)){
hes5.wilcox[i,2] <-  assign(paste(samples[i], ".wilcox", sep=""), wilcox.test(x=get(paste(samples[i],".B", sep=""))[,4], y=get(paste(samples[i],".T", sep=""))[,4], alternative="l", paired=TRUE))[3]
}
hist((-log2(hes5.wilcox[,2])), breaks=100, col="dark grey", xlab="Wilcox test log2( p-values )", main="Bisulphite sequencing HES5 promoter\nTumour-Normal pairs (n=39)")
abline(v=-log2(0.05), lty=2, col="blue")
text(x=-log2(0.05)+10, y=20, labels="0.05", col="blue")
```

![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6.png) 
___________

**Code and plots for Figure 2d-e (Prostate TCGA 450k array data):**


```r
#304 Tumour and 49 Benign HES5 core prom boxplot:
boxplot(tcga[which(tcga$V6 %in% core.hes5$SPOT_ID),10] ~ tcga[which(tcga$V6 %in% core.hes5$SPOT_ID),3], main="TCGA Prostate Samples\nHES5 methylation profiles")
```

![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-71.png) 

```r
#Wilcox tests for 49 Tumour-Benign pairs and histogram summary plot:
tn.samples <- levels(as.factor(tcga.tn.n[,11]))
tcga.tn.n[,17] <- paste(tcga.tn.n[,11], tcga.tn.n[,16], sep="_")
hes5.tcga.wilcox <- data.frame(tn.samples, rep(1, 49))
for(i in 1:49){
(hes5.tcga.wilcox[i,2] <- (wilcox.test(x=
tcga.tn.n[which(tcga.tn.n$V17==paste(tn.samples[i],"_Normal Tissue", sep="")),2], y=tcga.tn.n[which(tcga.tn.n$V17==paste(tn.samples[i],"_Primary Tumour", sep="")),2], alternative="l", paired=TRUE)[3]))
}

hist(-log2(hes5.tcga.wilcox[,2]), breaks=10, main="TCGA Prostate Tumour-Normal pairs\nHES5 methylation profiles (n=49)", xlab="Wilcox text -log2( p-value )")
abline(v=-log2(0.05), lty=2, col="blue")
text(x=-log2(0.05)+1, y=12, labels="0.05", col="blue")
```

![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-72.png) 
___________

**Code and plots for Figure 2f:**


```r
#extract Benign and Tumour avg meth:
full.data.t <- full.data[which(full.data[,4]=="T"),]
full.data.b <- full.data[which(full.data[,4]=="B"),]

#filter to remove pairs with missing data
full.data.t.core <- full.data.t[which(full.data.t[,3] %in% full.data.b[,3]),]
full.data.b.core <- full.data.b[which(full.data.b[,3] %in% full.data.t[,3]),]

#generate data for BS-seq ROC curve:

bsx.bn.100 <- rep(0, 100)
for(i in 1:100){
    (bsx.bn.100[i] <- length(which(full.data.b.core[,2]>i)))
}

bsx.t.100 <- rep(0, 100)
for(i in 1:100){
    (bsx.t.100[i] <- length(which(full.data.t.core[,2]>i)))
}



#ROCplot overlay BSx and TCGA::
require(pracma)
```

```
## Loading required package: pracma
```

```
## Warning: package 'pracma' was built under R version 3.1.1
```

```r
#ROC curves for Bisulphite sequening of the HES5 promoter (n=39 Tum.-Norm. pairs) and TCGA 450k array data for the HES5 promoter (n=49 Tum.-Norm. pairs):
plot(bsx.bn.100/39, bsx.t.100/39, type="l", lty=1, xlab="specificity", ylab="sensitivity", main="HES5 promoter methylation")
lines(tcga.bn.100/49, tcga.t.100/49, lty=2)
points(x=0.076, y=0.84, col="red", pch=19)

legend(x=0.37, y=0.8, legend=c(paste("BS-seq AUC = ",round(-1*trapz(bsx.bn.100/39,bsx.t.100/39), digits=3)),paste("TCGA AUC = ",round(-1*trapz(tcga.bn.100/49, tcga.t.100/49), digits=3)), paste("Pos.Pred.Val. = ", round(0.84/(0.84+0.076), digits=3))), lty=c(1,2,0), cex=0.8)
points(x=0.41, y=0.7, col="red", pch=19)
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8.png) 

```r
PPV<- round(0.85/(0.84+0.076), digits=3)
PPV
```

```
## [1] 0.928
```

```r
AUC.BS.seq <- round(-1*trapz(bsx.bn.100/39,bsx.t.100/39), digits=3)
AUC.BS.seq
```

```
## [1] 0.937
```

```r
AUC.TCGA.450k <- round(-1*trapz(tcga.bn.100/49, tcga.t.100/49), digits=3)
AUC.TCGA.450k
```

```
## [1] 0.905
```

___________
___________
**Figure 3**
-------------------------

**Code and plots for Figure 3a**

```r
load("cell-line-marmalaid.RData")

#library(marmalaid)
#library(minfi)
#data(annotation_v1.1)
#data(meth450)
#prostate <- annotation[which(annotation$TISSUE=="Prostate"),] 
#id <- as.character(prostate[,1])
#probes  <- row.names(meth450)
#beta <- getbeta(id,probes)
#status <- as.character(prostate[,8])

#get HES5 probe profiles:
#hes5.probes <- PriceAnno[which(PriceAnno$Closest_TSS_gene_name=="HES5"),22]
#gstp1.probes <- PriceAnno[which(PriceAnno$Closest_TSS_gene_name=="GSTP1"),22]

#beta.lines <- beta[,1:8]
#prostate.lines <- prostate[1:8,]

#re-order to group cell lines:
#beta.lines <- beta.lines[,c(1,2,3,7, 4, 5, 6, 8)]
#colnames(beta.lines) <- prostate.lines[,7]

#plot:
boxplot(beta.lines[hes5.probes[2:11],], las=2, ylab="Methylation proportion (Beta)", main="HES5 promoter")
```

![plot of chunk unnamed-chunk-9](figure/unnamed-chunk-91.png) 

```r
boxplot(beta.lines[gstp1.probes[2:11],], las=2, ylab="Methylation proportion (Beta)", main="GSTP1 promoter")
```

![plot of chunk unnamed-chunk-9](figure/unnamed-chunk-92.png) 
___________

**Code and plots for Figure 3b-c**



```r
load("lncap-aza-Gx.Rdata")

#library(GEOquery)
#load data:
#gse25346 <- getGEO(GEO="GSE25346", file="GSE25346_family.soft")

#load platform annots:
#gse25346.gpl <- GPLList(gse25346)
#GPL4133 <- Table(dataTable(gse25346.gpl$GPL4133))

#get data to single object:
#gse25346.gsm <- GSMList(gse25346)
#GSM623533 <- dataTable(gse25346.gsm$GSM623533)
#GSM623534 <- dataTable(gse25346.gsm$GSM623534)
#GSM623535 <- dataTable(gse25346.gsm$GSM623535)
#GSM623536 <- dataTable(gse25346.gsm$GSM623536)
#gse25346.df <- cbind(Table(GSM623533)[,2],Table(GSM623534)[,2],Table(GSM623535)[,2],Table(GSM623536)[,2]) 

#Boxplot genes on log2 scale:
#data downloaded as log10 ratios (LNCaP 24 or 48h aza/dmso - 2-color agilent arrays:)

boxplot(cbind(c(log2(10^(gse25346.df[which(GPL4133$GENE_SYMBOL=="HES5"),1:2]))), c(log2(10^(gse25346.df[which(GPL4133$GENE_SYMBOL=="HES5"),3:4])))), names=c("+aza-dC-24h", "+aza-dC-48h"), ylab="LNCaP log2( aza-dC / DMSO )", main="HES5 expression", ylim=c(-2.5, 2.5))
 abline(h=0, lty=2, col="light blue")
```

![plot of chunk unnamed-chunk-10](figure/unnamed-chunk-101.png) 

```r
boxplot(cbind(c(log2(10^(gse25346.df[which(GPL4133$GENE_SYMBOL=="HES6"),1:2]))), c(log2(10^(gse25346.df[which(GPL4133$GENE_SYMBOL=="HES6"),3:4])))), names=c("+aza-dC-24h", "+aza-dC-48h"), ylab="LNCaP log2( aza-dC / DMSO )", main="HES6 expression", ylim=c(-2.5, 2.5))
abline(h=0, lty=2, col="light blue")
```

![plot of chunk unnamed-chunk-10](figure/unnamed-chunk-102.png) 
___________

**Code and plots for Figure 3d-e**

```r
load("GSE3325-Varambally.RData")

###load geo.soft from pwd:
#gse3325 <- getGEO(GEO="GSE3325", file="GSE3325_family.soft")

###get platform info:
#gse3325.gpl <- GPLList(gse3325)
#GPL570 <- Table(dataTable(gse3325.gpl$GPL570))

###get data:
#gse3325.gsm <- GSMList(gse3325)

###make data.frame for all arrays:
#gse3325.df <- cbind(Table(gse3325.gsm[[1]])[,2], Table(gse3325.gsm[[2]])[,2], Table(gse3325.gsm[[3]])[,2], Table(gse3325.gsm[[4]])[,2], Table(gse3325.gsm[[5]])[,2], Table(gse3325.gsm[[6]])[,2], Table(gse3325.gsm[[7]])[,2], Table(gse3325.gsm[[8]])[,2], Table(gse3325.gsm[[9]])[,2], Table(gse3325.gsm[[10]])[,2], Table(gse3325.gsm[[11]])[,2], Table(gse3325.gsm[[12]])[,2], Table(gse3325.gsm[[13]])[,2], Table(gse3325.gsm[[14]])[,2], Table(gse3325.gsm[[15]])[,2], Table(gse3325.gsm[[16]])[,2], Table(gse3325.gsm[[17]])[,2], Table(gse3325.gsm[[18]])[,2], Table(gse3325.gsm[[19]])[,2])

###get sample names into usable format:
#samplesNames <- c(Meta(gse3325.gsm[[1]])$title, Meta(gse3325.gsm[[2]])$title, Meta(gse3325.gsm[[3]])$title, Meta(gse3325.gsm[[4]])$title, Meta(gse3325.gsm[[5]])$title, Meta(gse3325.gsm[[6]])$title, Meta(gse3325.gsm[[7]])$title, Meta(gse3325.gsm[[8]])$title, Meta(gse3325.gsm[[9]])$title, Meta(gse3325.gsm[[10]])$title, Meta(gse3325.gsm[[11]])$title, Meta(gse3325.gsm[[12]])$title, Meta(gse3325.gsm[[13]])$title, Meta(gse3325.gsm[[14]])$title, Meta(gse3325.gsm[[15]])$title, Meta(gse3325.gsm[[16]])$title, Meta(gse3325.gsm[[17]])$title, Meta(gse3325.gsm[[18]])$title, Meta(gse3325.gsm[[19]])$title)

###annotate data.frame:
#rownames(gse3325.df) <- GPL570[,11]
#colnames(gse3325.df) <- c(samplesNames)

###BN vs T boxplot
par(mar=c(5,6,4,2))
GENE="HES5"
boxplot(cbind(log2(as.numeric(gse3325.df[which(GPL570[,11]==GENE),1:6])), log2(as.numeric(gse3325.df[which(GPL570[,11]==GENE),7:12]))), names=c("Benign", "Primary Tumour"),ylab=paste(GENE, as.character(GPL570[which(GPL570[,11]==GENE),1]),"\nlog2(expression)"), main=paste(GENE, " expression\n(Varambally, et al 2005)"))
```

![plot of chunk unnamed-chunk-11](figure/unnamed-chunk-111.png) 

```r
par(mar=c(5,6,4,2))
GENE="HES6"
boxplot(cbind(log2(as.numeric(gse3325.df[which(GPL570[,11]==GENE),1:6])), log2(as.numeric(gse3325.df[which(GPL570[,11]==GENE),7:12]))), names=c("Benign", "Primary Tumour"),ylab=paste(GENE, as.character(GPL570[which(GPL570[,11]==GENE),1]),"\nlog2(expression)"), main=paste(GENE, " expression\n(Varambally, et al 2005)"))
```

![plot of chunk unnamed-chunk-11](figure/unnamed-chunk-112.png) 

```r
###Alternative for multi probe genes:
#par(mar=c(5,6,4,2))
#PROBE=1
#GENE="AMACR"
#boxplot(cbind(log2(as.numeric(gse3325.df[which(GPL570[,11]==GENE),1:6][PROBE,])), log2(as.numeric(gse3325.df[which(GPL570[,11]==GENE),7:12][PROBE,]))), names=c("Benign", "Primary Tumour"),ylab=paste(GENE, as.character(GPL570[which(GPL570[,11]==GENE),1][PROBE]),"\nlog2(expression)"), main=paste(GENE, " expression\n(Varambally, et al 2005)"))
```

___________

**Code and plots for Figure 3f-i**

```r
load("camcap-allin.RData")

#for(i in 1:length(core)){camcap.indexes[i] <- print(grep(core[i], phenoData(camcap.tumours.frozen)[[1]]))}

#camcap.hes5.samples <- (phenoData(camcap.tumours.frozen)[[1]][camcap.indexes])

#camcap.hes5.gx <- (assayData(camcap.tumours.frozen)[[1]][,camcap.indexes])

#camcap.hes5.samples.n <- (phenoData(camcap.normals.frozen)[[1]][camcap.indexes.n])

#camcap.hes5.gx.n <- (assayData(camcap.normals.frozen)[[1]][,camcap.indexes.n])

#for(i in 2:43){index.gx.for.meth[i] <- print(which(names(camcap.hes5.gx.core.ann)[i] == (hes5.meth.anno.T$sample)))}

#alt soln for bn:
#for(i in 2:43){index.gx.for.meth.B[i-1] <- print(which(names(camcap.hes5.gx.core.ann.n)[i] == (hes5.meth.anno.B$sample)))}

#camcap.hes5.gx.core.ann[index.gx.for.meth,]

#tumour.gx <- (assayData(camcap.tumours.frozen)[[1]])
#normal.gx <- (assayData(camcap.normals.frozen)[[1]])


#create object for scatter plots:
#HES5 <-  tumour.gx[rownames(camcap.hes5.gx.core.ann[which(camcap.hes5.gx.core.ann[,1]=="HES5"),]),]

#HES6 <-  tumour.gx[rownames(camcap.hes5.gx.core.ann[which(camcap.hes5.gx.core.ann[,1]=="HES6"),]),]

#HES1 <-  tumour.gx[rownames(camcap.hes5.gx.core.ann[which(camcap.hes5.gx.core.ann[,1]=="HES1"),]),]

#ERG <-  tumour.gx[rownames(camcap.hes5.gx.core.ann[which(camcap.hes5.gx.core.ann[,1]=="ERG"),]),][1,]

#scatter.genes <- rbind(HES6, HES5, ERG, HES6, ERG, HES1, HES1, HES6)

#rownames(scatter.genes) <- c("HES6", "HES5", "ERG", "HES6", "ERG","HES1","HES1", "HES6")

#rm(tumour.gx)
#rm(normal.gx)
#rm(camcap.normals.frozen)
#rm(camcap.tumours.frozen)


for(i in c(2,4,6,8)){
par(mar=c(6,6,4,2))
plot(scatter.genes[i-1,],scatter.genes[i,], xlab=rownames(scatter.genes)[i-1], ylab=rownames(scatter.genes)[i], pch=19, cex=0.7, xlim=c(min(scatter.genes[i-1,]), max(scatter.genes[i-1,])), ylim=c(min(scatter.genes[1,]), max(scatter.genes[1,])))
}
```

![plot of chunk unnamed-chunk-12](figure/unnamed-chunk-121.png) ![plot of chunk unnamed-chunk-12](figure/unnamed-chunk-122.png) ![plot of chunk unnamed-chunk-12](figure/unnamed-chunk-123.png) ![plot of chunk unnamed-chunk-12](figure/unnamed-chunk-124.png) 

```r
for(i in c(2,4,6,8)){print("Gene_1   Gene_2   correlation")
    print(paste(rownames(scatter.genes)[i-1],rownames(scatter.genes)[i], round(cor(scatter.genes[i-1,], scatter.genes[i,], method="pearson"),digits=3), sep="     "))
}
```

```
## [1] "Gene_1   Gene_2   correlation"
## [1] "HES6     HES5     0.032"
## [1] "Gene_1   Gene_2   correlation"
## [1] "ERG     HES6     -0.281"
## [1] "Gene_1   Gene_2   correlation"
## [1] "ERG     HES1     0.647"
## [1] "Gene_1   Gene_2   correlation"
## [1] "HES1     HES6     -0.244"
```
___________
___________
**Figure 4**
-------------------------

**Code and plots for Figure 4a-d:**


```r
load("ARG-workspace2.RData")

tc.gene.list <- c("TMPRSS2", "HES1", "Hes6", "KLK3")

for(i in 1:4){
par(mfrow=c(1,2), mar=c(6,4,1,0))
GENE=tc.gene.list[i]
probe=1

plot(x=VC.t[1:6], y=VCaP.centered[which(VCaP.centered$Genes==GENE),1:6][probe,],
     col="grey", cex=0.4, pch=1, lwd=2, ylab=paste(GENE, rownames(VCaP.centered[which(VCaP.centered$Genes==GENE),1:6][probe,])), ylim=c(-2,2), type="b",
     xlim=c(0,24))
points(x=VC.t[7:16], y=VCaP.centered[which(VCaP.centered$Genes==GENE),7:16][probe,],
       col="orange", cex=0.4, pch=1, lwd=2, type="b")           
points(x=VC.t[17:26], y=VCaP.centered[which(VCaP.centered$Genes==GENE),17:26][probe,],
       col="orange", cex=0.4, pch=1, lwd=2, type="b")           
points(x=VC.t[27:50], y=VCaP.centered[which(VCaP.centered$Genes==GENE),27:50][probe,],
       col="red", cex=0.4, pch=1, lwd=2, type="b")           
points(x=VC.t[51:70], y=VCaP.centered[which(VCaP.centered$Genes==GENE),51:70][probe,],
       col="red", cex=0.4, pch=1, lwd=2, type="b")           

points(x=LN.t[1:6], y=LNCaP.centered[which(LNCaP.all[,2]==GENE),1:6][probe,], 
     col="grey", cex=0.4, pch=1, lwd=2, type="b")
points(x=LN.t[7:11], y=LNCaP.centered[which(LNCaP.all[,2]==GENE),7:11][probe,], 
      col="light blue", cex=0.4, pch=1, lwd=2, type="b")
points(x=LN.t[12:38], y=LNCaP.centered[which(LNCaP.all[,2]==GENE),12:38][probe,], 
             col="blue", cex=0.4, pch=1, lwd=2, type="b")
points(x=LN.t[39:46], y=LNCaP.centered[which(LNCaP.all[,2]==GENE),39:46][probe,], 
       col="blue", cex=0.4, pch=1, lwd=2, type="b")
      
par(mar=c(6,4,1,1))
barplot(c(mean(as.numeric(LNCaP.V1[which(LNCaP.V1$Gene_symbol==GENE),3:7][probe,]))-5, 
          mean(as.numeric(LNCaP.R1[which(LNCaP.R1$Gene_symbol==GENE),3:29][probe,]))-5, 
          mean(as.numeric(LNCaP.R2[which(LNCaP.R2$Gene_symbol==GENE),3:10][probe,]))-5, 
          mean(as.numeric(VCaP.V1[which(VCaP.V1$Gene_symbol==GENE),3:12][probe,]))-5, 
          mean(as.numeric(VCaP.V2[which(VCaP.V2$Gene_symbol==GENE),3:12][probe,]))-5, 
          mean(as.numeric(VCaP.R1[which(VCaP.R1$Gene_symbol==GENE),3:26][probe,]))-5, 
          mean(as.numeric(VCaP.R2[which(VCaP.R2$Gene_symbol==GENE),3:22][probe,]))-5), 
        space=c(0.3,0.3,0,0.3,0,0.3,0),ylim=c(-2,12), col=c("light blue","blue", "blue" ,"orange", "orange", "red", "red"),
        names.arg=c("LNCaP.V", "LNCaP.R1", "LNCaP.R2", "VCaP.V1", "VCaP.V2", "VCaP.R1", "VCaP.R2"),
        las=2)

}
```

![plot of chunk unnamed-chunk-13](figure/unnamed-chunk-131.png) ![plot of chunk unnamed-chunk-13](figure/unnamed-chunk-132.png) ![plot of chunk unnamed-chunk-13](figure/unnamed-chunk-133.png) ![plot of chunk unnamed-chunk-13](figure/unnamed-chunk-134.png) 

___________

**Code and plots for Figure 4e-f:**


```r
#changepoint barplot:
LN.cp <- c(2,0,12,4)
VC.cp <- c(2,5,7,14)
N.cp <- c("TMPRSS2", "HES1", "HES6", "KLK3")
barplot(VC.cp, names.arg=N.cp, las=2, ylab="Expression change point (hrs)", col="red")
```

![plot of chunk unnamed-chunk-14](figure/unnamed-chunk-141.png) 

```r
barplot(LN.cp, names.arg=N.cp, las=2, ylab="Expression change point (hrs)", col="blue")
```

![plot of chunk unnamed-chunk-14](figure/unnamed-chunk-142.png) 
___________

**Code and plots for Figure 4g-l:**


```r
#xyplots for cor:

xy.tc.list <- c("TMPRSS2", "HES1", "Hes6", "HES1", "Hes6", "TMPRSS2")

for(i in c(2,4,5)){
par(mfrow=c(1,2))
xGENE=xy.tc.list[i-1]
yGENE=xy.tc.list[i]
xprobe=1
yprobe=1

plot(x=as.numeric(VCaP.centered[which(VCaP.centered$Genes==xGENE),1:6][xprobe,]), y=as.numeric(VCaP.centered[which(VCaP.centered$Genes==yGENE),1:6][yprobe,]), col="grey", cex=c(1.1^(VC.t[1:6])/3), pch=1, lwd=2, xlab=paste(xGENE, rownames(VCaP.centered[which(VCaP.centered$Genes==xGENE),1:6][xprobe,])), ylab=paste(yGENE, rownames(VCaP.centered[which(VCaP.centered$Genes==yGENE),1:6][yprobe,])), ylim=c(-2,2), xlim=c(-2,2), type="b")

points(x=as.numeric(VCaP.centered[which(VCaP.centered$Genes==xGENE),7:16][xprobe,]), y=as.numeric(VCaP.centered[which(VCaP.centered$Genes==yGENE),7:16][yprobe,]), col="orange", cex=c(1.1^(VC.t[7:16])/3), pch=1, lwd=2, type="b")           

points(x=as.numeric(VCaP.centered[which(VCaP.centered$Genes==xGENE),17:26][xprobe,]), y=as.numeric(VCaP.centered[which(VCaP.centered$Genes==yGENE),17:26][yprobe,]), col="orange", cex=c(1.1^(VC.t[17:26])/3), pch=1, lwd=2, type="b")           

points(x=as.numeric(VCaP.centered[which(VCaP.centered$Genes==xGENE),27:50][xprobe,]), y=as.numeric(VCaP.centered[which(VCaP.centered$Genes==yGENE),27:50][yprobe,]), col="red", cex=c(1.1^(VC.t[27:50])/3), pch=1, lwd=2, type="b")           

points(x=as.numeric(VCaP.centered[which(VCaP.centered$Genes==xGENE),51:70][xprobe,]), y=as.numeric(VCaP.centered[which(VCaP.centered$Genes==yGENE),51:70][yprobe,]), col="red", cex=c(1.1^(VC.t[51:70])/3), pch=1, lwd=2, type="b")           

plot(x=as.numeric(LNCaP.centered[which(LNCaP.all[,2]==xGENE),1:6][xprobe,]), y=as.numeric(LNCaP.centered[which(LNCaP.all[,2]==yGENE),1:6][yprobe,]), xlab=paste(xGENE, rownames(LNCaP.centered[which(LNCaP.centered$Gene_symbol==xGENE),1:6][xprobe,])), ylab=paste(yGENE, rownames(LNCaP.centered[which(LNCaP.centered$Gene_symbol==yGENE),1:6][yprobe,])), ylim=c(-2,2), type="b", xlim=c(-2,2), col="grey", cex=c(1.1^(VC.t[1:6])/3), pch=1, lwd=2)

points(x=as.numeric(LNCaP.centered[which(LNCaP.all[,2]==xGENE),7:11][xprobe,]), y=as.numeric(LNCaP.centered[which(LNCaP.all[,2]==yGENE),7:11][yprobe,]), 
      col="light blue", cex=c(1.1^(VC.t[7:11])/3), pch=1, lwd=2, type="b")

points(x=as.numeric(LNCaP.centered[which(LNCaP.all[,2]==xGENE),12:38][xprobe,]), y=as.numeric(LNCaP.centered[which(LNCaP.all[,2]==yGENE),12:38][yprobe,]), 
             col="blue", cex=c(1.1^(VC.t[12:38])/3), pch=1, lwd=2, type="b")

points(x=as.numeric(LNCaP.centered[which(LNCaP.all[,2]==xGENE),39:46][xprobe,]), y=as.numeric(LNCaP.centered[which(LNCaP.all[,2]==yGENE),39:46][yprobe,]), 
       col="blue", cex=c(1.1^(VC.t[39:46])/3), pch=1, lwd=2, type="b")
      
}
```

![plot of chunk unnamed-chunk-15](figure/unnamed-chunk-151.png) ![plot of chunk unnamed-chunk-15](figure/unnamed-chunk-152.png) ![plot of chunk unnamed-chunk-15](figure/unnamed-chunk-153.png) 

```r
#scale legends:
plot(x=rep(0:24,4), y=(rep(1:4, each=25)), cex=c(1.1^(0:24)/3), pch=1, lwd=2, col=c(rep("yellow", 25), rep("orange", 25), rep("red", 25), rep("dark red", 25)), ylim=c(0.5, 4.5))
plot(x=rep(0:24,3), y=(rep(1:3, each=25)), cex=c(1.1^(0:24)/3), pch=1, lwd=2, col=c(rep("light blue", 25), rep("blue", 25), rep("dark blue", 25)), ylim=c(0.5, 3.5))
```

![plot of chunk unnamed-chunk-15](figure/unnamed-chunk-154.png) 
___________


___________

___________

**Supplementary Figures**
====================================


**Supplementary Figure 1**
-------------------------

**Supplementary Figure 2b-2e: candidate region methylation boxplots**


```r
library(vioplot)
```

```
## Loading required package: sm
## Package 'sm', version 2.2-5.4: type help(sm) for summary information
```

```r
par(las=2)
vioplot(B006[,1], B006[,2], B006[,3], B006[,4], B006[,5], B006[,6], B006[,7], col="light blue", names=c(colnames(B006)))
```

```
## Warning: weights overwritten by binning
## Warning: weights overwritten by binning
## Warning: weights overwritten by binning
## Warning: weights overwritten by binning
## Warning: weights overwritten by binning
## Warning: weights overwritten by binning
## Warning: weights overwritten by binning
```

```r
title("B006 Methylation profile (B-value)")
```

![plot of chunk unnamed-chunk-16](figure/unnamed-chunk-161.png) 

```r
par(las=2)
vioplot(B007[,1], B007[,2], B007[,3], B007[,4], B007[,5], B007[,6], B007[,7], B007[,8], col="light blue", names=c(colnames(B007)))
```

```
## Warning: weights overwritten by binning
## Warning: weights overwritten by binning
## Warning: weights overwritten by binning
## Warning: weights overwritten by binning
## Warning: weights overwritten by binning
## Warning: weights overwritten by binning
## Warning: weights overwritten by binning
## Warning: weights overwritten by binning
```

```r
title("B007 Methylation profile (B-value)")
```

![plot of chunk unnamed-chunk-16](figure/unnamed-chunk-162.png) 

```r
par(las=2)
vioplot(B008[,1], B008[,2], B008[,3], B008[,4], col="light blue", names=c(colnames(B008)))
```

```
## Warning: weights overwritten by binning
## Warning: weights overwritten by binning
## Warning: weights overwritten by binning
## Warning: weights overwritten by binning
```

```r
title("B008 Methylation profile (B-value)")
```

![plot of chunk unnamed-chunk-16](figure/unnamed-chunk-163.png) 

```r
par(las=2)
vioplot(B054[,1], B054[,2], col="light blue", names=c("T1", "Benign"))
```

```
## Warning: weights overwritten by binning
## Warning: weights overwritten by binning
```

```r
title("B054 Methylation profile (B-value)")
```

![plot of chunk unnamed-chunk-16](figure/unnamed-chunk-164.png) 

**Additional to Supplementary Figure 2b-2e: Spearman rank correlation of B-values and binned methylation profiles**


```r
#B-value correlations
for(i in 1:6){print(cor(B006[,i], B006[,7], method="spear"))}
```

```
## [1] 0.9329
## [1] 0.9363
## [1] 0.9424
## [1] 0.9499
## [1] 0.9681
## [1] 0.9615
```

```r
for(i in 1:7){print(cor(B007[,i], B007[,8], method="spear"))}
```

```
## [1] 0.9584
## [1] 0.9542
## [1] 0.9644
## [1] 0.9621
## [1] 0.9646
## [1] 0.9722
## [1] 0.9681
```

```r
for(i in 1:2){print(cor(B008[,i], B008[,4], method="spear"))}
```

```
## [1] 0.9604
## [1] 0.9748
```

```r
cor(B054[,1], B054[,2], method="spear")
```

```
## [1] 0.9705
```

```r
#profile correlations:
for(i in 1:6){print(cor(as.numeric(hist(B006[,7], plot=FALSE)[[2]]), as.numeric(hist(B006[,i], plot=FALSE)[[2]])))}
```

```
## [1] 0.9666
## [1] 0.9958
## [1] 0.9615
## [1] 0.9842
## [1] 0.9971
## [1] 0.9815
```

```r
for(i in 1:7){print(cor(as.numeric(hist(B007[,8], plot=FALSE)[[2]]), as.numeric(hist(B007[,i], plot=FALSE)[[2]])))}
```

```
## [1] 0.9875
## [1] 0.9597
## [1] 0.9968
## [1] 0.9961
## [1] 0.987
## [1] 0.9919
## [1] 0.9838
```

```r
for(i in 1:3){print(cor(as.numeric(hist(B008[,4], plot=FALSE)[[2]]), as.numeric(hist(B008[,i], plot=FALSE)[[2]])))}
```

```
## [1] 0.9556
## [1] 0.9897
## [1] 0.9723
```

```r
cor(as.numeric(hist(B054[,2], plot=FALSE)[[2]]), as.numeric(hist(B054[,1], plot=FALSE)[[2]]))
```

```
## [1] 0.9422
```

**Supplementary Figure 2f-2m: candidate region methylation boxplots**


```r
gene.list <- rbind(c("HES5","chr1", 2461684, 2460184,"NM_001010926"), 
c("APC", "chr5", 112073556, 112181936, "NM_001127510"), 
c("RARB", "chr3",25469834, 25639422 , "NM_000965" ), 
c("TACC2", "chr10",123922941, 124014060, "NM_001291878"), 
c("ITGB2", "chr21", 46340965, 46305868, "NM_000211" ),
c("DGKZ", "chr11", 46383145, 46402104, "NM_001105540"), 
c("C5orf4", "chr5" , 154230213, 154198052, "NM_032385"), 
c("mir10b", "chr2",177015031, 177015140, "NR_029609"))

for(i in 1:dim(gene.list)[1]){

GENE <- gene.list[i,1]
CHR <- gene.list[i,2]
START <- as.numeric(gene.list[i,3])-200
END <- as.numeric(gene.list[i,3])+200

par(mfrow=c(1,4), mar=c(6,3.75,2,0), mgp=c(2.5,1,0))
boxplot(cbind(B006[which(PriceAnno$MAPINFO.1>=START & PriceAnno$MAPINFO.1<=END),],BBlood[which(PriceAnno$MAPINFO.1>=START & PriceAnno$MAPINFO.1<=END),1]), las=2, main="Case-006", names=c(colnames(B006[which(PriceAnno$MAPINFO.1>=START & PriceAnno$MAPINFO.1<=END),]), "Blood"), ylab=paste(GENE,"promoter methylation (B-value)"), border=col6ii, ylim=c(0,1))

par(mar=c(6,3.5,2,0))
boxplot(cbind(B007[which(PriceAnno$MAPINFO.1>=START & PriceAnno$MAPINFO.1<=END),],BBlood[which(PriceAnno$MAPINFO.1>=START & PriceAnno$MAPINFO.1<=END),2]), las=2, main="Case-007", names=c(colnames(B007[which(PriceAnno$MAPINFO.1>=START & PriceAnno$MAPINFO.1<=END),]), "Blood"), border=col7ii, ylim=c(0,1))

par(mar=c(6,3.5,2,2))
boxplot(cbind(B008[which(PriceAnno$MAPINFO.1>=START & PriceAnno$MAPINFO.1<=END),],BBlood[which(PriceAnno$MAPINFO.1>=START & PriceAnno$MAPINFO.1<=END),3]), las=2, main="Case-008", names=c(colnames(B008[which(PriceAnno$MAPINFO.1>=START & PriceAnno$MAPINFO.1<=END),]), "Blood"), border=col8ii, ylim=c(0,1))

par(mar=c(6,1.5,2,7.5))
boxplot(B054[which(PriceAnno$MAPINFO.1>=START & PriceAnno$MAPINFO.1<=END),], las=2, names=c("T1", "Benign"), main="Case-054", border=col54ii, ylim=c(0,1))

}
```

![plot of chunk unnamed-chunk-18](figure/unnamed-chunk-181.png) ![plot of chunk unnamed-chunk-18](figure/unnamed-chunk-182.png) ![plot of chunk unnamed-chunk-18](figure/unnamed-chunk-183.png) ![plot of chunk unnamed-chunk-18](figure/unnamed-chunk-184.png) ![plot of chunk unnamed-chunk-18](figure/unnamed-chunk-185.png) ![plot of chunk unnamed-chunk-18](figure/unnamed-chunk-186.png) ![plot of chunk unnamed-chunk-18](figure/unnamed-chunk-187.png) ![plot of chunk unnamed-chunk-18](figure/unnamed-chunk-188.png) 

Alternative plots showing Tumour / Benign (to highlight tumour specificity of methylation changes):


```r
gene.list <- rbind(c("HES5","chr1", 2461684, 2460184,"NM_001010926"), 
                   c("APC", "chr5", 112073556, 112181936, "NM_001127510"), 
                   c("RARB", "chr3",25469834, 25639422 , "NM_000965" ), 
                   c("TACC2", "chr10",123922941, 124014060, "NM_001291878"), 
                   c("ITGB2", "chr21", 46340965, 46305868, "NM_000211" ),
                   c("DGKZ", "chr11", 46383145, 46402104, "NM_001105540"), 
                   c("C5orf4", "chr5" , 154230213, 154198052, "NM_032385"), 
                   c("mir10b", "chr2",177015031, 177015140, "NR_029609"))

for(i in 1:dim(gene.list)[1]){
  
  GENE <- gene.list[i,1]
  CHR <- gene.list[i,2]
  START <- as.numeric(gene.list[i,3])-200
  END <- as.numeric(gene.list[i,3])+200
  
  par(mfrow=c(1,4), mar=c(6,3.75,4,0), mgp=c(2.5,1,0))
  boxplot(B006[which(PriceAnno$MAPINFO.1>=START & PriceAnno$MAPINFO.1<=END),1:6]/B006[which(PriceAnno$MAPINFO.1>=START & PriceAnno$MAPINFO.1<=END),7], las=2, main=paste(GENE, "\nCase-006"), ylab="Tumour/Benign", border=col6ii, ylim=c(0,25), names=c(colnames(B006[,1:6])))
  
  par(mar=c(6,3.5,4,0))
  boxplot(B007[which(PriceAnno$MAPINFO.1>=START & PriceAnno$MAPINFO.1<=END),1:7]/B007[which(PriceAnno$MAPINFO.1>=START & PriceAnno$MAPINFO.1<=END),8], las=2, main=paste(GENE, "\nCase-007"), border=col7ii, ylim=c(0,25), names=c(colnames(B007[,1:7])))
  
  par(mar=c(6,3.5,4,2))
  boxplot(B008[which(PriceAnno$MAPINFO.1>=START & PriceAnno$MAPINFO.1<=END),1:3]/B008[which(PriceAnno$MAPINFO.1>=START & PriceAnno$MAPINFO.1<=END),4], las=2, main=paste(GENE, "\nCase-008"), border=col8ii, ylim=c(0,25), names=c(colnames(B008[,1:3])))
  
  par(mar=c(6,1.5,4,7.5))
  boxplot(B054[which(PriceAnno$MAPINFO.1>=START & PriceAnno$MAPINFO.1<=END),1]/B054[which(PriceAnno$MAPINFO.1>=START & PriceAnno$MAPINFO.1<=END),2], las=2, main=paste(GENE, "\nCase-054"), border=col54ii, ylim=c(0,25), names=c(colnames(B054[,1])))
}
```

![plot of chunk unnamed-chunk-19](figure/unnamed-chunk-191.png) ![plot of chunk unnamed-chunk-19](figure/unnamed-chunk-192.png) ![plot of chunk unnamed-chunk-19](figure/unnamed-chunk-193.png) ![plot of chunk unnamed-chunk-19](figure/unnamed-chunk-194.png) ![plot of chunk unnamed-chunk-19](figure/unnamed-chunk-195.png) ![plot of chunk unnamed-chunk-19](figure/unnamed-chunk-196.png) ![plot of chunk unnamed-chunk-19](figure/unnamed-chunk-197.png) ![plot of chunk unnamed-chunk-19](figure/unnamed-chunk-198.png) 
___________

**Supplementary Figure 2**
-------------------------



```r
gene.list <- rbind(c("HES5","chr1", 2461684, 2460184,"NM_001010926"), 
c("APC", "chr5", 112073556, 112181936, "NM_001127510"), 
c("RARB", "chr3",25469834, 25639422 , "NM_000965" ), 
c("TACC2", "chr10",123922941, 124014060, "NM_001291878"), 
c("ITGB2", "chr21", 46340965, 46305868, "NM_000211" ),
c("DGKZ", "chr11", 46383145, 46402104, "NM_001105540"), 
c("C5orf4", "chr5" , 154230213, 154198052, "NM_032385"), 
c("mir10b", "chr2",177015031, 177015140, "NR_029609"))

for(i in 1:dim(gene.list)[1]){
GENE.bucket <- gene.list[i,]
WIN <- 10000

GENE <- GENE.bucket[1]
CHR <- GENE.bucket[2]
START <- as.numeric(GENE.bucket[3])
END <- as.numeric(GENE.bucket[4])
NAME <- GENE.bucket[5]
#STRAND <- GENE.bucket[6]
#START1 <- if(STRAND=="+"){print(START)}else{END}
#END1 <- if(STRAND=="-"){print(START)}else{END}

#tidy-up with GENE.PROBE.ROWS object
GENE.PROBE.ROWS <- which(seqnames(hm450)==CHR & start(hm450)>START-WIN & end(hm450)<START+WIN)

par(mar=c(7,6,4,2), mgp=c(5,0.7,0), mfrow=c(1,1))
plot(x=(PriceAnno[GENE.PROBE.ROWS,4]), y=log(B006[GENE.PROBE.ROWS,1]/B006[GENE.PROBE.ROWS,7]), type="b", col="red", main=paste(GENE, " locus"),
 cex=0.3, xlab=paste("Genomic coordinates of ",GENE," probes"), ylab="Methylation change\nlog(Tumour/benign)", las=2, ylim=c(-2,4))
   points(x=(PriceAnno[GENE.PROBE.ROWS,4]), y=log(B006[GENE.PROBE.ROWS,2]/B006[GENE.PROBE.ROWS,7]), type="b", col="red", cex=0.3)
   points(x=(PriceAnno[GENE.PROBE.ROWS,4]), y=log(B006[GENE.PROBE.ROWS,3]/B006[GENE.PROBE.ROWS,7]), type="b", col="red", cex=0.3)
   points(x=(PriceAnno[GENE.PROBE.ROWS,4]), y=log(B006[GENE.PROBE.ROWS,4]/B006[GENE.PROBE.ROWS,7]), type="b", col="light blue", cex=0.3)
   points(x=(PriceAnno[GENE.PROBE.ROWS,4]), y=log(B006[GENE.PROBE.ROWS,5]/B006[GENE.PROBE.ROWS,7]), type="b", col="light green", cex=0.3)
   points(x=(PriceAnno[GENE.PROBE.ROWS,4]), y=log(B006[GENE.PROBE.ROWS,6]/B006[GENE.PROBE.ROWS,7]), type="b", col="orange", cex=0.3)
   points(x=(PriceAnno[GENE.PROBE.ROWS,4]), y=log(B006[GENE.PROBE.ROWS,7]/B006[GENE.PROBE.ROWS,7]), type="b", col="black", cex=0.3)

lines(x=c(START,END), y=c(-1.5, -1.5))
text(x=START+500, y=-1.9, labels=NAME)
legend(x=START+200, y=4, legend=c(colnames(B006)), lty=1, col=col6, cex=0.6)

}
```

![plot of chunk unnamed-chunk-20](figure/unnamed-chunk-201.png) ![plot of chunk unnamed-chunk-20](figure/unnamed-chunk-202.png) ![plot of chunk unnamed-chunk-20](figure/unnamed-chunk-203.png) ![plot of chunk unnamed-chunk-20](figure/unnamed-chunk-204.png) ![plot of chunk unnamed-chunk-20](figure/unnamed-chunk-205.png) ![plot of chunk unnamed-chunk-20](figure/unnamed-chunk-206.png) ![plot of chunk unnamed-chunk-20](figure/unnamed-chunk-207.png) ![plot of chunk unnamed-chunk-20](figure/unnamed-chunk-208.png) 

__________

**Supplementary Figure 3**
-------------------------


```r
#CpG density accross HES5 prom. amplicon:
hist(as.numeric(rownames(table(N702_N507.CpG))), xlab="HES5 promotor amplicon position (bp)", ylab="CpGs per 50bp window", main="CpG density across HES5 promotor amplicon")
```

![plot of chunk unnamed-chunk-21](figure/unnamed-chunk-211.png) 

```r
plot(x=rownames(table(N707_N506.CpG)), y=((rowSums(table(N707_N506.CpG)))), type="h", xlab="postition of CpG in HES5 promoter", ylab="Sequencing coverage", main="Benign TB09.1008")
```

![plot of chunk unnamed-chunk-21](figure/unnamed-chunk-212.png) 

```r
#plot coverage at each CpG:
plot(x=rownames(table(N702_N507.CpG)), y=((rowSums(table(N702_N507.CpG)))), type="h", xlab="postition of CpG in HES5 promoter", ylab="Sequencing coverage", main="Tumour TB09.1008")
```

![plot of chunk unnamed-chunk-21](figure/unnamed-chunk-213.png) 


__________

**Supplementary Figure 4**
-------------------------


```r
load("hes5Validation2.RData")

for(i in 1:length(samples)){
par(mfrow=c(2,1), mar=c(4,5,1.3,2), mgp=c(2,0.5,0))
plot(get(paste(samples[i],".B", sep=""))[,4], type='h', ylim=c(0,50), xlab="HES5 prom. CpGs", ylab="% methylation", main=paste("Benign ",samples[i], " bisulphite seq.", sep=""))
plot(get(paste(samples[i],".T", sep=""))[,4], type='h', ylim=c(0,50), xlab="HES5 prom. CpGs", ylab="% methylation", main=paste("Tumour ",samples[i], " bisulphite seq.", sep=""))
}
```

![plot of chunk unnamed-chunk-22](figure/unnamed-chunk-221.png) ![plot of chunk unnamed-chunk-22](figure/unnamed-chunk-222.png) ![plot of chunk unnamed-chunk-22](figure/unnamed-chunk-223.png) ![plot of chunk unnamed-chunk-22](figure/unnamed-chunk-224.png) ![plot of chunk unnamed-chunk-22](figure/unnamed-chunk-225.png) ![plot of chunk unnamed-chunk-22](figure/unnamed-chunk-226.png) ![plot of chunk unnamed-chunk-22](figure/unnamed-chunk-227.png) ![plot of chunk unnamed-chunk-22](figure/unnamed-chunk-228.png) ![plot of chunk unnamed-chunk-22](figure/unnamed-chunk-229.png) ![plot of chunk unnamed-chunk-22](figure/unnamed-chunk-2210.png) ![plot of chunk unnamed-chunk-22](figure/unnamed-chunk-2211.png) ![plot of chunk unnamed-chunk-22](figure/unnamed-chunk-2212.png) ![plot of chunk unnamed-chunk-22](figure/unnamed-chunk-2213.png) ![plot of chunk unnamed-chunk-22](figure/unnamed-chunk-2214.png) ![plot of chunk unnamed-chunk-22](figure/unnamed-chunk-2215.png) ![plot of chunk unnamed-chunk-22](figure/unnamed-chunk-2216.png) ![plot of chunk unnamed-chunk-22](figure/unnamed-chunk-2217.png) ![plot of chunk unnamed-chunk-22](figure/unnamed-chunk-2218.png) ![plot of chunk unnamed-chunk-22](figure/unnamed-chunk-2219.png) ![plot of chunk unnamed-chunk-22](figure/unnamed-chunk-2220.png) ![plot of chunk unnamed-chunk-22](figure/unnamed-chunk-2221.png) ![plot of chunk unnamed-chunk-22](figure/unnamed-chunk-2222.png) ![plot of chunk unnamed-chunk-22](figure/unnamed-chunk-2223.png) ![plot of chunk unnamed-chunk-22](figure/unnamed-chunk-2224.png) ![plot of chunk unnamed-chunk-22](figure/unnamed-chunk-2225.png) ![plot of chunk unnamed-chunk-22](figure/unnamed-chunk-2226.png) ![plot of chunk unnamed-chunk-22](figure/unnamed-chunk-2227.png) ![plot of chunk unnamed-chunk-22](figure/unnamed-chunk-2228.png) ![plot of chunk unnamed-chunk-22](figure/unnamed-chunk-2229.png) ![plot of chunk unnamed-chunk-22](figure/unnamed-chunk-2230.png) ![plot of chunk unnamed-chunk-22](figure/unnamed-chunk-2231.png) ![plot of chunk unnamed-chunk-22](figure/unnamed-chunk-2232.png) ![plot of chunk unnamed-chunk-22](figure/unnamed-chunk-2233.png) ![plot of chunk unnamed-chunk-22](figure/unnamed-chunk-2234.png) ![plot of chunk unnamed-chunk-22](figure/unnamed-chunk-2235.png) ![plot of chunk unnamed-chunk-22](figure/unnamed-chunk-2236.png) ![plot of chunk unnamed-chunk-22](figure/unnamed-chunk-2237.png) ![plot of chunk unnamed-chunk-22](figure/unnamed-chunk-2238.png) ![plot of chunk unnamed-chunk-22](figure/unnamed-chunk-2239.png) 

__________

**Supplementary Figure 5**
-------------------------
**Supplementary Figure 5c-5d: expression xyplots of HES5-vs-HES6 and HES1-vs-HES6 with point size scaled by ERG expression (to allow 3d-visulaization of ERG, HES1, HES6 relationship)**


```r
load("camcap-allin.RData")

plot(HES1, HES6, xlab=paste(rownames(camcap.hes5.gx.core.ann[which(camcap.hes5.gx.core.ann[,1]=="HES1"),]), "\nHES1"), ylab=paste(rownames(camcap.hes5.gx.core.ann[which(camcap.hes5.gx.core.ann[,1]=="HES6"),]), "\nHES6"), pch=19, cex=ERG-5.7, main=paste("Point size scaled by\nERG:", rownames(camcap.hes5.gx.core.ann[which(camcap.hes5.gx.core.ann[,1]=="ERG"),])[1]))
```

![plot of chunk unnamed-chunk-23](figure/unnamed-chunk-231.png) 

```r
plot(HES5, HES6, xlab=paste(rownames(camcap.hes5.gx.core.ann[which(camcap.hes5.gx.core.ann[,1]=="HES5"),]), "\nHES5"), ylab=paste(rownames(camcap.hes5.gx.core.ann[which(camcap.hes5.gx.core.ann[,1]=="HES6"),]), "\nHES6"), pch=19, cex=ERG-5.7, main=paste("Point size scaled by\nERG:", rownames(camcap.hes5.gx.core.ann[which(camcap.hes5.gx.core.ann[,1]=="ERG"),])[1]))
```

![plot of chunk unnamed-chunk-23](figure/unnamed-chunk-232.png) 

```r
plot(y=ERG, x=rep(1,length(ERG)),cex=ERG-5.7, pch=19, ylab="ERG expression")
```

![plot of chunk unnamed-chunk-23](figure/unnamed-chunk-233.png) 

__________

**Supplementary Figure 7**
-------------------------

**Supplementary Figure 7c-h**


```r
#xyplots for cor:

xy.tc.list <- c("HES1","ERG","HES5", "ERG", "HES5", "TMPRSS2")

for(i in c(2,4,5)){
par(mfrow=c(1,2))
xGENE=xy.tc.list[i-1]
yGENE=xy.tc.list[i]
xprobe=1
yprobe=1

plot(x=as.numeric(VCaP.centered[which(VCaP.centered$Genes==xGENE),1:6][xprobe,]), y=as.numeric(VCaP.centered[which(VCaP.centered$Genes==yGENE),1:6][yprobe,]), col="grey", cex=c(1.1^(VC.t[1:6])/3), pch=1, lwd=2, xlab=paste(xGENE, rownames(VCaP.centered[which(VCaP.centered$Genes==xGENE),1:6][xprobe,])), ylab=paste(yGENE, rownames(VCaP.centered[which(VCaP.centered$Genes==yGENE),1:6][yprobe,])), ylim=c(-2,2), xlim=c(-2,2), type="b")

points(x=as.numeric(VCaP.centered[which(VCaP.centered$Genes==xGENE),7:16][xprobe,]), y=as.numeric(VCaP.centered[which(VCaP.centered$Genes==yGENE),7:16][yprobe,]), col="orange", cex=c(1.1^(VC.t[7:16])/3), pch=1, lwd=2, type="b")           

points(x=as.numeric(VCaP.centered[which(VCaP.centered$Genes==xGENE),17:26][xprobe,]), y=as.numeric(VCaP.centered[which(VCaP.centered$Genes==yGENE),17:26][yprobe,]), col="orange", cex=c(1.1^(VC.t[17:26])/3), pch=1, lwd=2, type="b")           

points(x=as.numeric(VCaP.centered[which(VCaP.centered$Genes==xGENE),27:50][xprobe,]), y=as.numeric(VCaP.centered[which(VCaP.centered$Genes==yGENE),27:50][yprobe,]), col="red", cex=c(1.1^(VC.t[27:50])/3), pch=1, lwd=2, type="b")           

points(x=as.numeric(VCaP.centered[which(VCaP.centered$Genes==xGENE),51:70][xprobe,]), y=as.numeric(VCaP.centered[which(VCaP.centered$Genes==yGENE),51:70][yprobe,]), col="red", cex=c(1.1^(VC.t[51:70])/3), pch=1, lwd=2, type="b")           

plot(x=as.numeric(LNCaP.centered[which(LNCaP.all[,2]==xGENE),1:6][xprobe,]), y=as.numeric(LNCaP.centered[which(LNCaP.all[,2]==yGENE),1:6][yprobe,]), xlab=paste(xGENE, rownames(LNCaP.centered[which(LNCaP.centered$Gene_symbol==xGENE),1:6][xprobe,])), ylab=paste(yGENE, rownames(LNCaP.centered[which(LNCaP.centered$Gene_symbol==yGENE),1:6][yprobe,])), ylim=c(-2,2), type="b", xlim=c(-2,2), col="grey", cex=c(1.1^(VC.t[1:6])/3), pch=1, lwd=2)

points(x=as.numeric(LNCaP.centered[which(LNCaP.all[,2]==xGENE),7:11][xprobe,]), y=as.numeric(LNCaP.centered[which(LNCaP.all[,2]==yGENE),7:11][yprobe,]), 
      col="light blue", cex=c(1.1^(VC.t[7:11])/3), pch=1, lwd=2, type="b")

points(x=as.numeric(LNCaP.centered[which(LNCaP.all[,2]==xGENE),12:38][xprobe,]), y=as.numeric(LNCaP.centered[which(LNCaP.all[,2]==yGENE),12:38][yprobe,]), 
             col="blue", cex=c(1.1^(VC.t[12:38])/3), pch=1, lwd=2, type="b")

points(x=as.numeric(LNCaP.centered[which(LNCaP.all[,2]==xGENE),39:46][xprobe,]), y=as.numeric(LNCaP.centered[which(LNCaP.all[,2]==yGENE),39:46][yprobe,]), 
       col="blue", cex=c(1.1^(VC.t[39:46])/3), pch=1, lwd=2, type="b")
      
}
```

![plot of chunk unnamed-chunk-24](figure/unnamed-chunk-241.png) ![plot of chunk unnamed-chunk-24](figure/unnamed-chunk-242.png) ![plot of chunk unnamed-chunk-24](figure/unnamed-chunk-243.png) 
