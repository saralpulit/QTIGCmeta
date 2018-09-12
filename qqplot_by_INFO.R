#
# This script generates a QQ plot after reading a file that contains z scores or chi-sq values
# Usage (p-values): R CMD BATCH -input.pvals -PVAL -output.pdf qqplot.R
#
#

## READS input options
rm(list=ls())

x <- 0 

repeat {
        x <- x+1
        if (commandArgs()[x] == "-CL") {
        input <- commandArgs()[x+1]; input <- substr(input, 2, nchar(input))
        stat_type <- commandArgs()[x+2]; stat_type <- substr(stat_type, 2, nchar(stat_type))
        output <- commandArgs()[x+3]; output <- substr(output, 2, nchar(output))
        break
        }
        if (x == length(commandArgs())) {
                print("remember the -CL command!")
                break}
        }

rm(x)

print(input)

#input=commandArgs()[7]
#input=substr(input,2,nchar(input))
#
#output=commandArgs()[9]
#output=substr(output,2,nchar(output))

## Plot function ##
plotQQ <- function(z,color,cex){
p <- 2*pnorm(-abs(z))
p <- sort(p)
expected <- c(1:length(p))
lobs <- -(log10(p))
lexp <- -(log10(expected / (length(expected)+1)))

# plots all points with p < 1e-3
p_sig = subset(p,p<0.001)

points(lexp[1:length(p_sig)], lobs[1:length(p_sig)], pch=23, cex=.3, col=color, bg=color)

# samples 2500 points from p > 1e-3
n=2501
i<- c(length(p)- c(0,round(log(2:(n-1))/log(n)*length(p))),1)
lobs_bottom=subset(lobs[i],lobs[i] <= 3)
lexp_bottom=lexp[i[1:length(lobs_bottom)]]

if(length(lobs_bottom) != length(lexp_bottom)) {
	lobs_bottom <- lobs_bottom[1:length(lexp_bottom)]
}

print(length(lobs_bottom))
print(length(lexp_bottom))

points(lexp_bottom, lobs_bottom, pch=23, cex=cex, col=color, bg=color)
}


## Reads data
S <- read.table(input,header=F)
#pvals=S[,1]

pvals_lo1=subset(S, ( S[,2] >= 0.0 & S[,2] <= 0.20 ))
pvals_lo2=subset(S, ( S[,2] > 0.20 & S[,2] <= 0.40 ))
pvals_lo3=subset(S, ( S[,2] > 0.40 & S[,2] <= 0.60 ))
pvals_lo4=subset(S, ( S[,2] > 0.60 & S[,2] <= 0.80 ))
pvals_lo5=subset(S, ( S[,2] > 0.80 & S[,2] <= 1.00 ))
pvals_lo6=subset(S, ( S[,2] > 1.00 ))

z=qnorm(S[,1]/2)
z_lo1=qnorm(pvals_lo1[,1]/2)
z_lo2=qnorm(pvals_lo2[,1]/2)
z_lo3=qnorm(pvals_lo3[,1]/2)
z_lo4=qnorm(pvals_lo4[,1]/2)
z_lo5=qnorm(pvals_lo5[,1]/2)
z_lo6=qnorm(pvals_lo6[,1]/2)

print(z_lo6)

## calculates lambda
lambda = round(median(z^2)/qchisq(0.5,df=1),3)

## Plots axes and null distribution
pdf(output, width=6, height=6)
plot(c(0,7), c(0,7), col="red", lwd=3, type="l", xlab="Expected Distribution (-log10 of P value)", ylab="Observed Distribution (-log10 of P value)", xlim=c(0,7), ylim=c(0,7), las=1, xaxs="i", yaxs="i", bty="l",main=c(substitute(paste("QQ plot: ",lambda," = ", lam),list(lam = lambda)),expression()))

## plots data

plotQQ(z,"black",0.4);
plotQQ(z_lo1,"olivedrab1",0.3);
plotQQ(z_lo2,"orange",0.3);
plotQQ(z_lo3,"lightskyblue",0.3);
plotQQ(z_lo4,"purple",0.3);
plotQQ(z_lo5,"darkgreen",0.3);
plotQQ(z_lo6,"darkblue",0.3);

## provides legend
legend(.25,7,legend=c("Expected (null)","Observed",
paste("0.0 < INFO < 0.2 [",length(z_lo1),"]"),
paste("0.2 < INFO < 0.4 [",length(z_lo2),"]"),
paste("0.4 < INFO < 0.6 [",length(z_lo3),"]"),
paste("0.6 < INFO < 0.8 [",length(z_lo4),"]"),
paste("0.8 < INFO < 1.0 [",length(z_lo5),"]"),
paste("1.0 > INFO [",length(z_lo6),"]")),
pch=c((vector("numeric",5)+1)*23), cex=c((vector("numeric",5)+0.8)), pt.bg=c("red","black","olivedrab1","orange","lightskyblue", "purple", "darkgreen", "darkblue"))
rm(z)
dev.off()
