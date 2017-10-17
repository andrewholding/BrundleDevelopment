######################
#
# Settings
#
######################

setwd("/Volumes/FlyPeaks/FlyPeaks")

######################
#
# Functions 
#
######################

source('package/brundle.R')

######################
#
# Main Code
#
######################

filename<-"Rdata/040_SLX-14438_dba_human_ER_CTCF.rda"
load(filename)

pdf("plots/045_SLX-14438_CTCF_Fulvestrant.pdf")
dba.plotVenn(dbaControl,dbaControl$masks$Fulvestrant, label1="Rep1", label2="Rep2", label3="Rep3", main="CTCF Peaks Fulvestrant")
    dev.off()
pdf("plots/045_SLX-14438_CTCF_Control.pdf")
    dba.plotVenn(dbaControl,dbaControl$masks$none, label1="Rep1", label2="Rep2", label3="Rep3", main="CTCF Peaks Control")
dev.off()

pdf("plots/045_SLX-14438_ER_Fulvestrant.pdf")
    dba.plotVenn(dbaExperiment,dbaExperiment$masks$Fulvestrant, label1="Rep1", label2="Rep2", label3="Rep3", main="ER Peaks Fulvestrant")
dev.off()
pdf("plots/045_SLX-14438_ER_control.pdf")
    dba.plotVenn(dbaExperiment,dbaExperiment$masks$none,label1="Rep1", label2="Rep2", label3="Rep3", main="ER Peaks Control")
dev.off()
    
filename<-"Rdata/026_SLX-8047_dba_drosophila.rda"
load(filename)

pdf("plots/045_SLX-8047_H2av_Fulvestrant.pdf")
    dba.plotVenn(dba_dm,dba_dm$masks$Fulvestrant, label1="Rep1", label2="Rep2", label3="Rep3",label4="Rep4", main="H2av Peaks Fulvestrant")
dev.off()
pdf("plots/045_SLX-8047_H2av_control.pdf")
    dba.plotVenn(dba_dm,dba_dm$masks$none,label1="Rep1", label2="Rep2", label3="Rep3",label4="Rep4", main="H2av Peaks Control")
dev.off()

filename<-"Rdata/026_SLX-8047_dba_human.rda"
load(filename)

pdf("plots/045_SLX-8047_ER_Fulvestrant.pdf")
    dba.plotVenn(dba,dba$masks$Fulvestrant, label1="Rep1", label2="Rep2", label3="Rep3", label4="Rep4",main="ER Peaks Fulvestrant")
dev.off()
pdf("plots/045_SLX-8047_ER_Control.pdf")
    dba.plotVenn(dba,dba$masks$none,label1="Rep1", label2="Rep2", label3="Rep3",label4="Rep4", main="ER Peaks Control")
dev.off()
    
filename<-"Rdata/014_SLX-12998_dba_mouse.rda"
load(filename)

pdf("plots/045_SLX-12998_mmER_Fulvestrant.pdf")
dba.plotVenn(dba,dba$masks$Fulvestrant, label1="Rep1", label2="Rep2", label3="Rep3", label4="Rep4",main="Mouse ER Peaks Fulvestrant")
dev.off()
pdf("plots/045_SLX-12998_mmER_Control.pdf")
    dba.plotVenn(dba,dba$masks$none,label1="Rep1", label2="Rep2", label3="Rep3",label4="Rep4", main="Mouse ER Peaks Control")
dev.off()

filename<-"Rdata/015_SLX-12998_dba_human.rda"
load(filename)

pdf("plots/045_SLX-12998_hsER_Fulvestrant.pdf")
    dba.plotVenn(dba,dba$masks$Fulvestrant, label1="Rep1", label2="Rep2", label3="Rep3", label4="Rep4",main="ER Peaks Fulvestrant")
dev.off()
pdf("plots/045_SLX-12998_hsER_Control.pdf")
    dba.plotVenn(dba,dba$masks$none,label1="Rep1", label2="Rep2", label3="Rep3",label4="Rep4", main="ER Peaks Control")
dev.off()