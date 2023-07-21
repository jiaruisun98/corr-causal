library(openxlsx)
library(ggplot2)
library(dplyr)
library(patchwork)

###############################
#This code uses simulation results in Simulation.R and Simulation_regression.R to produce Figure 1, 2 and S1.
#The simulation results are first organized in the xlsx form, see `Simulation-result.xlsx`
#plot_block produces Figure 1
#plot_off produces Figure 2
#plot_lasso produces Figure S1


#####Plot-block-correlation#####

Result_block<-read.xlsx("Simulation-result.xlsx",sheet=1)
Result_block<-as.data.frame(Result_block)
Result_block$rho<-as.factor(Result_block$rho)
Result_block$dependency<-as.factor(Result_block$dependency)
levels(Result_block$dependency)<-c("rho[time]~'='~0.3", "rho[time]~'='~0")

plot_block<-ggplot(data=Result_block,aes(x=n_p,y=value,color=Method_plot,shape=rho,linetype=rho,group=Method))+
  geom_point(size=3)+
  geom_line(size=1)+
  theme_bw()+
  facet_grid(dependency~FOM,labeller=label_parsed)+
  scale_color_manual(values = c('blue','red'),name="Method")+
  scale_linetype_manual(values = c("solid","dashed","dotted"))+
  scale_y_continuous(breaks=c(0.0,0.05,0.2,0.4,0.6,0.8))+
  geom_hline(yintercept = 0.05,linetype="dotted",color="black",show.legend =TRUE)+ 
  labs(x="Sample size and dimension",y="Value")+
  theme(axis.title = element_text(size = 16))+
  theme(strip.text.x = element_text(size = 14, colour = "black"))+
  theme(strip.text.y = element_text(size = 14, colour = "black"))+
  theme(legend.position = "bottom")+
  theme(legend.text = element_text(size=12))+
  theme(legend.title = element_text(size=12))
plot_block$labels$linetype<-expression(rho[paste(0," ")])
plot_block$labels$shape<-expression(rho[paste(0," ")])

plot_block

######Plot-off diagonal-correlation#####


Result_off<-read.xlsx("Simulation-result.xlsx",sheet=2)
Result_off<-as.data.frame(Result_off)
Result_off$rho<-as.factor(Result_off$rho)
Result_off$dependency<-as.factor(Result_off$dependency)
levels(Result_off$dependency)<-c("rho[time]~'='~0.3", "rho[time]~'='~0")



plot_off<-ggplot(data=Result_off,aes(x=n_p,y=value,color=Method_plot,shape=rho,linetype=rho,group=Method))+
  geom_point(size=3)+
  geom_line(size=1)+
  theme_bw()+
  facet_grid(dependency~FOM,labeller=label_parsed)+
  scale_color_manual(values = c('blue','red'),name="Method")+
  scale_linetype_manual(values = c("solid","dashed","dotted"))+
  scale_y_continuous(breaks=c(0.0,0.05,0.2,0.4,0.6,0.8))+
  geom_hline(yintercept = 0.05,linetype="dotted",color="black",show.legend =TRUE)+ 
  labs(x="Sample size and dimension",y="Value")+
  theme(axis.title = element_text(size = 16))+
  theme(strip.text.x = element_text(size = 14, colour = "black"))+
  theme(strip.text.y = element_text(size = 14, colour = "black"))+
  theme(legend.position = "bottom")+
  theme(legend.text = element_text(size=12))+
  theme(legend.title = element_text(size=12))
plot_off$labels$linetype<-expression(rho[paste(0," ")])
plot_off$labels$shape<-expression(rho[paste(0," ")])

plot_off


######Plot-block-regression#####



Result_lasso<-read.xlsx("Simulation-result.xlsx",sheet=3)
Result_lasso<-as.data.frame(Result_lasso)
Result_lasso$rho<-as.factor(Result_lasso$rho)
Result_lasso$dependency<-as.factor(Result_lasso$dependency)
levels(Result_lasso$dependency)<-c("rho[time]~'='~0.3", "rho[time]~'='~0")



plot_lasso<-ggplot(data=Result_lasso,aes(x=n_p,y=value,color=Method_plot,shape=rho,linetype=rho,group=Method))+
  geom_point(size=3)+
  geom_line(size=1)+
  theme_bw()+
  facet_grid(dependency~FOM,labeller=label_parsed)+
  scale_color_manual(values = c('red','blue'),name="Method")+
  scale_linetype_manual(values = c("solid","dashed","dotted"))+
  scale_y_continuous(breaks=c(0.0,0.05,0.2,0.4,0.6,0.8))+
  geom_hline(yintercept = 0.05,linetype="dotted",color="black",show.legend =TRUE)+ 
  labs(x="Sample size and dimension",y="Value")+
  theme(axis.title = element_text(size = 16))+
  theme(strip.text.x = element_text(size = 14, colour = "black"))+
  theme(strip.text.y = element_text(size = 14, colour = "black"))+
  theme(legend.position = "bottom")+
  theme(legend.text = element_text(size=12))+
  theme(legend.title = element_text(size=12))
plot_lasso$labels$linetype<-expression(rho[paste(0," ")])
plot_lasso$labels$shape<-expression(rho[paste(0," ")])

plot_lasso























