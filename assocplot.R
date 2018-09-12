# 
# This script generates a regional association plot
#
# Author: Paul de Bakker, debakker@broad.mit.edu
#
# Last update: 20 December 2010 (SLP)
#
# Usage : R CMD BATCH -input.txt -genes.txt -rs2880058 -48002 -100 -output.pdf assocplot.R
#
#
#options(echo=FALSE)

## READS input options
rm(list=ls())
input=commandArgs()[7]
input=substr(input,2,nchar(input))

input2=commandArgs()[8]
input2=substr(input2,2,nchar(input2))

index=commandArgs()[9]
index=substr(index,2,nchar(index))

n.eff=commandArgs()[10]
n.eff=substr(n.eff,2,nchar(n.eff))

range.y=commandArgs()[11]
range.y=substr(range.y,2,nchar(range.y))

locusname=commandArgs()[12]
locusname=substr(locusname,2,nchar(locusname))

output=commandArgs()[13]
output=substr(output,2,nchar(output))

#library(Cairo)

#make.fancy.locus.plot(index, genes, locusname, locus[1,1], locus, as.numeric(range.y), as.numeric(n.eff), display.r2=T, left.axis=T, right.axis=T, ld.hue=0)

make.fancy.locus.plot <- function(snp, genes, locusname, chr, locus, range.y, n.eff, size.pos = 1000000, min.pos = -1, display.r2 = T, left.axis = T, right.axis = T, ld.hue=0) {

hit <- locus[snp,]
#hit2 <- locus["rs11708996",]

#
# size of the region
#
#min.pos <- min(locus$POS) - 10000
#max.pos <- max(locus$POS) + 10000
#size.pos <- max.pos - min.pos
#center.pos <- min.pos + ( size.pos / 2 )
#center.100kb.pos <- round(center.pos / 100000) * 100000
#offset.100kb.pos <- round((size.pos/3) / 100000) * 100000

if ( min.pos == -1 ) {
  center.pos <- hit$POS
  min.pos <- hit$POS - ( size.pos / 2 )
  max.pos <- hit$POS + ( size.pos / 2 )
} 
else {
  center.pos <- min.pos + ( size.pos / 2 )
  max.pos <- min.pos + size.pos 
}

center.100kb.pos <- round(center.pos / 100000) * 100000
offset.100kb.pos <- round((size.pos/3) / 100000) * 100000

#
# range of y-axis
#
# this dedicates 33% of the yaxis to the genes, labels, recomb rate

if ( range.y < 0 ) { 
  range.y <- max( 6, ceiling( -(log10( hit$PVALUE )) / 10 ) * 10 * 1.10 )
}

### IF PVALUE IS ZERO, ARBITRARILY SET TO 1E-150 FOR VISUALIZATION

if ( hit$PVALUE == 0 ) {
  range.y <- max( 6, ceiling( -(log10( 1e-150 )) / 10 ) * 10 * 1.10 )
}


range.y <- range.y*1.5

offset <- range.y / 3 
big.range <- range.y + offset 

ystart.gene <- - offset
#ystart.recomb <- - offset + (big.range / 8)
ystart.recomb <- 0

#
# recombination rate 
#
recomb <- read.table(paste("/fg/software/Tagger/assocplot/genetic_map_chr", chr, "_b36.txt", sep=""), header=T)
keep.recomb <- subset(recomb, recomb[,1] > min.pos & recomb[,1] < max.pos)

# Default r^2 to zero if not known
locus$RSQR[is.na(locus$RSQR)] <- 0


#
# genes in the region
#
#genelist <- read.table("refseq_genes_032109.short.largest_footprint.txt", header=T)
#genes.in.locus <- subset(genelist, ( genelist$CHR == as.character(chr) & ( ( genelist$START > min.pos & genelist$START < max.pos ) | ( genelist$STOP > min.pos & genelist$STOP < max.pos) ) ) )
#print(genes.in.locus)

#
# genotyped markers
#
#markers.hip <- subset(locus, (row.names(locus) != snp & locus$POS > min.pos & locus$POS < max.pos & locus$PVALUE > 1e-3 & ! is.na(locus$RSQR) ))
#markers.medp <- subset(locus, (row.names(locus) != snp & locus$POS > min.pos & locus$POS < max.pos & locus$PVALUE < 1e-3 & locus$PVALUE > 1e-6 & ! is.na(locus$RSQR) ))
#markers.lowp <- subset(locus, (row.names(locus) != snp & locus$POS > min.pos & locus$POS < max.pos & locus$PVALUE < 1e-6 & locus$PVALUE > 5e-8 & ! is.na(locus$RSQR) ))
#markers.ld <- subset(locus, (row.names(locus) != snp & locus$POS > min.pos & locus$POS < max.pos & ! is.na(locus$RSQR) ))
#markers.nold <- subset(locus, (row.names(locus) != snp & locus$POS > min.pos & locus$POS < max.pos & is.na(locus$RSQR) ))

markers.nsig <- subset(locus, (row.names(locus) != snp & locus$POS > min.pos & locus$POS < max.pos & locus$PVALUE > 5e-8 & is.na(locus$FUNCTION) ))
markers.nsig.missense <- subset(locus, (row.names(locus) != snp & locus$POS > min.pos & locus$POS < max.pos & locus$PVALUE > 5e-8 & locus$FUNCTION == "missense" ))

markers.sigp <- subset(locus, (row.names(locus) != snp & locus$POS > min.pos & locus$POS < max.pos & locus$PVALUE <= 5e-8 & is.na(locus$FUNCTION) ))
markers.sigp.missense <- subset(locus, (row.names(locus) != snp & locus$POS > min.pos & locus$POS < max.pos & locus$PVALUE <= 5e-8 & locus$FUNCTION == "missense" ))


par(mar=c(4,4,3,4))

#
# start plot with recombination rate (in background)
#
plot(keep.recomb[,1], ystart.recomb + ( ( keep.recomb[,2] / 90 ) * ( 3 * big.range / 4 )), type="l", col="lightblue", lwd=1, xlim=c(min.pos, max.pos), ylim=c(-offset,range.y), xlab="", ylab="", main=locusname, axes=F)


#
# axes, titles and legends
#
mtext(paste("Chromosome", chr, "position (kb)", sep=" "), side=1, line=2.5)
axis(1, at=c(center.100kb.pos - offset.100kb.pos, center.100kb.pos, center.100kb.pos + offset.100kb.pos), labels=c((center.100kb.pos - offset.100kb.pos) / 1000, center.100kb.pos / 1000, (center.100kb.pos + offset.100kb.pos) / 1000), las=1) 

#
# left-hand axis
#
if ( left.axis == T ) {
  #axis(2, at=seq(0,range.y,2), labels=seq(0,range.y,2), las=1) 
  #per.tick <- round( (range.y/10), 0 ) * 2
  per.tick <- ceiling(range.y/5)
  axis(2, at=seq(0,range.y,per.tick), labels=seq(0,range.y,per.tick), las=1) 
  mtext(expression(paste(-log[10]," of P-value")), side=2, at=(range.y/2), line=2.5)
}

#
# right-hand axis
#
if ( right.axis == T ) {
  axis(4, at=c( ystart.recomb, ystart.recomb + (big.range / 4), ystart.recomb + ( 2 * big.range / 4), ystart.recomb + ( 3 * big.range / 4 ) ), labels=c("0","30","60","90"), las=1)
  #mtext("Recombination rate (cM/Mb)", side=4, at=(-offset+big.range/2), line=2)
  text(par("usr")[2], range.y, srt=-90, pos = 4, offset = 3, labels = "Recombination rate (cM/Mb)", xpd =TRUE)
}

#box()
#lines(c(min.pos, max.pos), c(0,0), lty="dotted", lwd=1, col="black")
lines(c(min.pos, max.pos), c(-(log10(5e-8)),-(log10(5e-8))), lty="dashed", lwd=1, col="black")


#
# plot the markers
#
points(markers.nsig$POS, -(log10(markers.nsig$PVALUE)), pch=23, cex=0.6, bg=hsv(ld.hue,markers.nsig$RSQR,1), col="grey20", lwd=0.2)
points(markers.nsig.missense$POS, -(log10(markers.nsig.missense$PVALUE)), pch=25, cex=0.6, bg=hsv(ld.hue,markers.nsig.missense$RSQR,1), col="darkblue", lwd=1.0)
points(markers.sigp$POS, -(log10(markers.sigp$PVALUE)), pch=23, cex=sqrt(markers.sigp$N_EFF/n.eff), bg=hsv(ld.hue,markers.sigp$RSQR,1), col="grey20", lwd=0.2)
points(markers.sigp.missense$POS, -(log10(markers.sigp.missense$PVALUE)), pch=25, cex=sqrt(markers.sigp.missense$N_EFF/n.eff), bg=hsv(ld.hue,markers.sigp.missense$RSQR,1), col="darkblue", lwd=1.0)

#
# this is the hit
#
#points(hit$POS, -(log10(hit$PVALUE)), pch=23, cex=2.5, bg="red", col=rgb(1-markers.hip$N_EFF/n.eff,1-markers.hip$N_EFF/n.eff,1-markers.hip$N_EFF/n.eff))
#points(hit$POS, -(log10(hit$PVALUE)), pch=23, cex=hit$N_EFF/n.eff, bg="red", col="black", lwd=0.3)

points(hit$POS, -(log10(hit$PVALUE)), pch=23, cex=1.5, bg=hsv(ld.hue,1,1), col="black", lwd=0.3)
text(hit$POS, -(log10(hit$PVALUE)), labels=c(row.names(hit)), pos=3, cex=0.6)


### IF THE PVALUE IS ZERO, ARBITRARILY SET TO 1E-150 FOR VISUALIZATION

if ( hit$PVALUE == 0 ) {
  points(hit$POS, -(log10( 1e-150 )), pch=23, cex=1.5, bg=hsv(ld.hue,1,1), col="black", lwd=0.3)
  text(hit$POS, -(log10( 1e-150 )), labels=c(row.names(hit)), pos=3, cex=0.6)
}

#points(hit2$POS, -(log10(hit2$PVALUE)), pch=23, cex=1.5, bg=hsv(ld.hue,1,1), col="black", lwd=0.3)
#text(hit2$POS, -(log10(hit2$PVALUE)), labels=c(row.names(hit2)), pos=3, cex=0.6)

#if ( -(log10(best.pval)) < range.y ) {
#	points(hit$POS, -(log10(best.pval)), pch=23, cex=2.5, bg="blue")
#	text(hit$POS, -(log10(best.pval)), labels=c(paste("P=",best.pval,sep="")), pos=4, offset=2)
#}
#else {
#	points(hit$POS, range.y, pch=23, cex=2.5, bg="blue")
#	text(hit$POS, range.y, labels=c(paste("P=",best.pval,sep="")), pos=4, offset=1)
#}
#

#
# r^2 legend
#
if ( display.r2 == T ) {
  step.y <- range.y / 20
  x.pos.r2 <- max.pos - ( size.pos * 0.005 )

  points(x.pos.r2, range.y - (step.y * 11), pch=22, cex=2.5, bg=hsv(ld.hue,0,1), col=hsv(ld.hue,0,1))
  points(x.pos.r2, range.y - (step.y * 10), pch=22, cex=2.5, bg=hsv(ld.hue,0.1,1), col=hsv(ld.hue,0.1,1))
  points(x.pos.r2, range.y - (step.y * 9), pch=22, cex=2.5, bg=hsv(ld.hue,0.2,1), col=hsv(ld.hue,0.2,1))
  points(x.pos.r2, range.y - (step.y * 8), pch=22, cex=2.5, bg=hsv(ld.hue,0.3,1), col=hsv(ld.hue,0.3,1))
  points(x.pos.r2, range.y - (step.y * 7), pch=22, cex=2.5, bg=hsv(ld.hue,0.4,1), col=hsv(ld.hue,0.4,1))
  points(x.pos.r2, range.y - (step.y * 6), pch=22, cex=2.5, bg=hsv(ld.hue,0.5,1), col=hsv(ld.hue,0.5,1))
  points(x.pos.r2, range.y - (step.y * 5), pch=22, cex=2.5, bg=hsv(ld.hue,0.6,1), col=hsv(ld.hue,0.6,1))
  points(x.pos.r2, range.y - (step.y * 4), pch=22, cex=2.5, bg=hsv(ld.hue,0.7,1), col=hsv(ld.hue,0.7,1))
  points(x.pos.r2, range.y - (step.y * 3), pch=22, cex=2.5, bg=hsv(ld.hue,0.8,1), col=hsv(ld.hue,0.8,1))
  points(x.pos.r2, range.y - (step.y * 2), pch=22, cex=2.5, bg=hsv(ld.hue,0.9,1), col=hsv(ld.hue,0.9,1))
  points(x.pos.r2, range.y - step.y, pch=22, cex=2.5, bg=hsv(ld.hue,1.0,1), col=hsv(ld.hue,1.0,1))

  text(x.pos.r2, range.y - (step.y * 2), labels=c(paste("0.8")), pos=1, cex=0.6)
  text(x.pos.r2, range.y - (step.y * 5), labels=c(paste("0.5")), pos=1, cex=0.6)
  text(x.pos.r2, range.y - (step.y * 11), labels=c(expression(paste(r^2))), pos=1)
}


###
### draw some genes
###
  track.ypos <- offset / 4 
  name.offset <- offset / 10

  for (i in 1:nrow(genes) ) { 
      if ( ( genes[i,]$START * 1000 > min.pos & genes[i,]$START * 1000 < max.pos ) || 
           ( genes[i,]$STOP * 1000 > min.pos & genes[i,]$STOP * 1000 < max.pos ) ) {
           y.pos <- -offset + ( genes[i,]$TRACK * track.ypos )
           arrows( max(genes[i,]$START * 1000, min.pos), y.pos, min(genes[i,]$STOP * 1000, max.pos), y.pos, length=0.025, code=genes[i,]$CODE, lwd=2, angle=45, lty="solid", col="darkgreen", lend="butt", ljoin="mitre")
           if ( genes[i,]$NAMEPOS * 1000 < min.pos ) {
               text( min.pos, y.pos + name.offset, labels=genes[i,]$NAME, cex=0.6, font=3, pos=4 )
           } else if ( genes[i,]$NAMEPOS * 1000 > max.pos ) {
               text( max.pos, y.pos + name.offset, labels=genes[i,]$NAME, cex=0.6, font=3, pos=4 )
           } else {
               text( genes[i,]$NAMEPOS * 1000, y.pos + name.offset, labels=genes[i,]$NAME, cex=0.6, font=3, pos=4 ) 
           }
      }
  }       

}

locus <- read.table(input, header=T, row.names=1)

genes <- read.table(input2, header=T)



pdf(output, width=6, height=4)

make.fancy.locus.plot(index, genes, locusname, locus[1,1], locus, as.numeric(range.y), as.numeric(n.eff), display.r2=T, left.axis=T, right.axis=T, ld.hue=0)

#for IGFBP3
#make.fancy.locus.plot(index, genes, locusname, locus[1,1], locus, as.numeric(range.y), as.numeric(n.eff), min.pos=45800000)

#for CDKN2C
#make.fancy.locus.plot(index, genes, locusname, locus[1,1], locus, as.numeric(range.y), as.numeric(n.eff), size.pos=1500000)

#for PLN
#make.fancy.locus.plot(index, genes, locusname, locus[1,1], locus, as.numeric(range.y), as.numeric(n.eff), min.pos=118400000)

#for SCN5A
#make.fancy.locus.plot(index, genes, locusname, locus[1,1], locus, as.numeric(range.y), as.numeric(n.eff), min.pos=38000000, size.pos=1500000)
#make.fancy.locus.plot(index, genes, locusname, locus[1,1], locus, as.numeric(range.y), as.numeric(n.eff), min.pos=38000000, size.pos=1500000, ld.hue=0.7)
#make.fancy.locus.plot(index, genes, locusname, locus[1,1], locus, as.numeric(range.y), as.numeric(n.eff), min.pos=38000000, size.pos=1500000, ld.hue=0.4)


dev.off()

quit()

