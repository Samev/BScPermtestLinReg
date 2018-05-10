#File used to plot results for permutation tests, as currently configured to be placed on the same folder as the result files

#For processing results
processRes <- function(dir=getwd(), sign, df, beta, n, reps, perms, tests, dataType){
  res <- matrix(nrow=length(df)*length(beta)*length(n), ncol=length(tests)+1+5)
  count <- 0
  for(i in 1:length(beta)){
    for (j in 1:length(df)){
      for (k in 1:length(n)){
        count <- count + 1
        dataParams <- c(n[k], beta[i], df[j])
        filename <- paste(dir, dataType, gsub("\\.", "", paste(dataParams, collapse="-")), "_n", n[k], "_reps", reps, "_perms", perms, "_", paste(tests, collapse="-"), ".rds", sep="")
        pres <- readRDS(filename)
        pres <- matrix(data=unlist(pres), nrow=6, ncol=500)
        pres[4,] <- abs(pres[4,])
        rejRate <- apply(pres[1:(length(tests)+1),] <= c(sign, sign, sign, qt(1-sign/2, n[k]-2)), 1, sum)/reps
        rejRate[4] <- 1-rejRate[4]
        kurt <- sum(pres[length(tests)+1+1,])/reps
        skew <- sum(pres[length(tests)+1+2,])/reps
        res[count, ] <- c(rejRate, kurt, skew, n[k], df[j], beta[i])
      }
    }
  }
  res <- data.frame(res)
  colnames(res) <- c("rSquaredTest", "tTest2min", "tTestabs", "normTest", "kurtosis", "skewness", "n", "v", "beta")
  return(res)
}

sign<-0.05 #Significance level
dataType <- "simpleExpSkew_XExp" #Specify data type
load("parameters.RData")

dir <- "./" # Set what directory the result files are found in
res <- processRes(dir, sign, paramVec, betaVec[betaVec!=0], nVec, reps, perms, testFuns, dataType)
res0 <- processRes(dir, sign, paramVec, 0, nVec, reps, perms, testFuns, dataType)

# Reference interval for Type1-error
RI <- data.frame(x=c(-Inf, Inf), ys=sign+sqrt(sign*(1-sign)/reps)*c(1.96,-1.96))



#################### Datatype-specific #########################################

nrOfBetas = 3 #Number of non-zero beta values

xaxis=c("skewness") # What is to be imported in ggplot (our x-axis), needs to exist in res and res0
xlabel = "Skevhet"  # The x-label of the plots

ys=c("tTest2min", "normTest") # Which data (y) that should be used, needs to exist in res and res0
shapesLab=c('"permTest"', '"t-test"')      # What these are to be called in the legend
# Note the format of the data on previous row, need to be exactly: '"name"'


################################################################################


#plot
lineSize <- 0.6
pointSize <- 2.5

require('ggplot2')
require('extrafont')
#font_import(pattern = "lmroman*") #to import the font (Computer Modern), make sure it is installed on the computer
#loadfonts()

for(i in nVec){
  plot1=ggplot()
  
  for(j in 1:length(ys)){
    plot1 = plot1+
      geom_line(data=res0[res0$n == i,], aes_string(x=xaxis, y=ys[j]), size=lineSize) +
      geom_point(data=res0[res0$n == i,], aes_string(x=xaxis, y=ys[j], shape=shapesLab[j]), size=pointSize)
  }
  
  plot1=plot1+
    geom_rect(data=RI, aes(xmin=x[1], xmax=x[2], ymin=ys[2], ymax=ys[1], fill="95% R.I, Typ1"), alpha=0.1)+
    scale_fill_manual(values=c("#999999"))+
    scale_color_manual(values=c("#222222","#222222"))+
    scale_shape_manual(values=(0:length(shapesLab)))+
    labs(fill="", shape='Testtyp: ',title=bquote(n == .(i)~", "~beta[1]==0))+
    ylim(0,0.10) +
    xlab(xlabel)+ylab("Förkastningsfrekvens")+
    guides(fill = guide_legend(order=1))+
    theme_bw()+
    theme(plot.title = element_text(hjust = 0.5),
          legend.position="bottom", text=element_text(size=10, family="LM Roman 10"),
          legend.spacing.x = unit(0.1, 'cm'),
          legend.spacing.y = unit(0.05, 'cm'),
          legend.key.height = unit(0.1,"cm"),
          legend.key.width = unit(0.5,"cm"),
          legend.box = "vertical",
          legend.box.spacing = unit(0.1,"cm"))
  #print(plot1)
  
  name <- paste(dataType, "_beta0_n", i, sep="")
  ggsave(plot=plot1, filename=paste(name, "pdf", sep="."), path=dir, width = 8, height = 9, units = "cm", family="LM Roman 10")
  
  
  plot2=ggplot()
  
  for(j in 1:length(ys)){
    plot2 = plot2+ 
      geom_line(data=res[res$n == i,], aes_string(x=xaxis, y=ys[j], linetype='factor(beta)', color='factor(beta)'), size=lineSize) +
      geom_point(data=res[res$n == i,], aes_string(x=xaxis, y=ys[j], shape=shapesLab[j],color='factor(beta)'), size=pointSize)
  }
  
  plot2 = plot2+
    scale_linetype_manual(values = (nrOfBetas:1)) +
    scale_shape_manual(values=(0:length(shapesLab)))+
    scale_color_manual(values=c("#aaaaaa","#777777","#222222"))+
    labs(color=bquote(beta[1]), linetype=bquote(beta[1]), shape="Testtyp :", title=bquote(n == .(i)~", "~beta[1]!=0))+
    ylim(0,1) +
    xlab(xlabel)+ylab("Förkastningsfrekvens")+
    guides(linetype = guide_legend(order=1), color = guide_legend(order=1))+
    theme_bw()+
    theme(plot.title = element_text(hjust = 0.5),
          legend.position="bottom", text=element_text(size=10, family="LM Roman 10"),
          legend.spacing.x = unit(0.1, 'cm'),
          legend.spacing.y = unit(0.05, 'cm'),
          legend.key.height = unit(0.1,"cm"),
          legend.key.width = unit(0.5,"cm"),
          legend.box = "vertical",
          legend.box.spacing = unit(0.1,"cm"))
  #  print(plot2)
  name <- paste(dataType, "_betaN0_n", i, sep="")
  ggsave(plot=plot2, filename=paste(name, "pdf", sep="."), path=dir, width = 8, height = 9, units = "cm", family="LM Roman 10")
}
