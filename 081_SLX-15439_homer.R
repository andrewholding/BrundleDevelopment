#Homer profile

ctcf<-read.table("SLX-15439/blacklist_filtered/Homer_CTCF_Profile.txt", sep="\t", header=TRUE)
er<-read.table("SLX-15439/blacklist_filtered/Homer_ER_Profile.txt", sep="\t", header=TRUE)


par(mfrow=c(2,2))
profile_CTCF<-as.matrix(ctcf[,c(1,2,5,8,11)])
plot(profile_CTCF,t="l",lwd=3,col="#00AFBB", xlab="Distance from Summit",ylab="Read Depth",main="CTCF Profile")
lines(profile_CTCF[,c(1,3)],lwd=3,col="#E7B800")
lines(profile_CTCF[,c(1,4)],lwd=3,col="#FC4E07")
lines(profile_CTCF[,c(1,5)],lwd=3,col="#FC00FE")
legend("topright",legend = c("PDX01", "PDX03", "PDX04", "PDX05"), col=c("#00AFBB","#E7B800", "#FC4E07","#FC00FE"),pch = c(19,19))

maxCTCF<-colMaxs(profile_CTCF)
minCTCF<-colMins(profile_CTCF)

profile_CTCF[,2:5]<-profile_CTCF[,2:5]-minCTCF[2:5]

colSums(profile_CTCF[,2:5]/maxCTCF[2:5])
normalised_CTCF<-t((t(profile_CTCF[,2:5])/maxCTCF[2:5]))
normalised_CTCF<-cbind(profile_CTCF[,1],normalised_CTCF)



plot(normalised_CTCF,lwd=3,t="l",col="#00AFBB", xlab="Distance from Summit",ylab="Read Depth",main="Normalised CTCF Profile")
lines(normalised_CTCF[,c(1,3)],lwd=3,col="#E7B800",)
lines(normalised_CTCF[,c(1,4)],lwd=3,col="#FC4E07")
lines(normalised_CTCF[,c(1,5)],lwd=3,col="#FC00FE")
legend("topright",legend = c("PDX01", "PDX03", "PDX04", "PDX05"), col=c("#00AFBB","#E7B800", "#FC4E07","#FC00FE"),pch = c(19,19))




profile_ER<-as.matrix(er[,c(1,2,5,8,11)])
plot(profile_ER,t="l",lwd=3,col="#00AFBB", xlab="Distance from Summit",ylab="Read Depth",main="ER Profile")
lines(profile_ER[,c(1,3)],lwd=3,col="#E7B800")
lines(profile_ER[,c(1,4)],lwd=3,col="#FC4E07")
lines(profile_ER[,c(1,5)],lwd=3,col="#FC00FE")
legend("topright",legend = c("PDX01", "PDX03", "PDX04", "PDX05"), col=c("#00AFBB","#E7B800", "#FC4E07","#FC00FE"),pch = c(19,19))



normalised_ER<-t((t(profile_ER[,2:5])/maxCTCF[2:5]))
normalised_ER<-cbind(profile_ER[,1],normalised_ER)


plot(normalised_ER,t="l",lwd=3,col="#00AFBB", xlab="Distance from Summit",ylab="Read Depth",main="Normalised ER Profile")
lines(normalised_ER[,c(1,3)],lwd=3,col="#E7B800")
lines(normalised_ER[,c(1,4)],lwd=3,col="#FC4E07")
lines(normalised_ER[,c(1,5)],lwd=3,col="#FC00FE")
legend("topright",legend = c("PDX01", "PDX03", "PDX04", "PDX05"), col=c("#00AFBB","#E7B800", "#FC4E07","#FC00FE"),pch = c(19,19))

