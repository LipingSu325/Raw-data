# GSE116256 PMID: 33648605
# GSE161901 PMID: 33227818
# GSE148218 PMID: 33218352
# GSE37642


if (T) {
  dir.create("data")
  dir.create("scripts")
  dir.create("results")
  dir.create("files")
  dir.create("PDFs")
  dir.create("00_origin_datas",recursive = T)
  
}
library(stringr)
library(openxlsx)
library(data.table)
library(reshape2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(clusterProfiler)
library(pheatmap)
library(ComplexHeatmap)
library(GSVA)
library(GSEABase)
library(fgsea)
library(corrplot)
library(colorspace)
library(survival)
library(survminer)
library(maftools)
library(vegan)
library(forcats)
library(ggpubr)
library(ggsci)
library(ggplot2)
library(rstatix)
library(ggstatsplot)
library(ggcor)
library(ggstance)
rm(list = ls())
options(stringsAsFactors = F)

# source('z:/projects/codes/mg_base.R')
my_mutiboxplot=function(dat,group,group_cols=ggsci::pal_aaas()(10),
                        test_method=c('t.test','wilcox.test','anova','kruskal.test')[2],
                        bw=T,xlab='',ylab='score',title='',size=10,angle = 45, hjust = 1,
                        legend.position='top',fill='group',notch=F){
  # dat=tcga.est[tcga.subtype.cli$Samples,]
  # group=tcga.subtype.cli$Cluster
  dat.bind=cbind(dat,Cluster=group)
  dat.bind=crbind2DataFrame(dat.bind)
  dat.melt=melt(dat.bind)
  #data=data[which(!is.na(data[,1])),]
  colnames(dat.melt)=c('Group','type','value')
  if(length(names(table(group)))>2){
    test_method='kruskal.test'
  }else{
    test_method='wilcox.test'
  }
  p=dat.melt %>%
    ggplot(aes(x=type, y=value,fill=Group)) +
    geom_boxplot(notch = notch) +  
    scale_fill_manual(values =group_cols)+   #
    ggpubr::stat_compare_means(aes(group=Group), label = "p.signif", method = test_method)+
    labs(x="", y = ylab, fill =fill,title =title) +
    #theme_light()+
    # theme_bw()+
    #theme_classic()
    theme(legend.position = legend.position,                 #
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(size = size,angle = angle, hjust = hjust)) # 
  return(p)
}
my_riskplot=function(cli_dat,cols=c("red","blue"),xlab='Samples',
                     a.ylab="Risk score",b.labs="Survival time(year)",cutoff=0,labs=c('A','B')){
  #cli_dat=tcga.risktype.cli
  cli.dat.order=cli_dat[order(cli_dat$Riskscore),c('OS.time','Status','Riskscore','Risktype')]
  fp_dat=data.frame(Samples=1:nrow(cli_dat),cli.dat.order)
  p1=ggplot(fp_dat,aes(x=Samples,y=Riskscore))+geom_point(aes(color=Risktype))+
    scale_colour_manual(values =cols)+
    theme_bw()+labs(x=xlab,y=a.ylab)+
    geom_hline(yintercept=cutoff,colour="black", linetype="dotted",size=0.8)+
    geom_vline(xintercept=sum(fp_dat$Risktype=="Low"),colour="black", linetype="dotted",size=0.8)
  p1
  p2=ggplot(fp_dat,aes(x=Samples,y=OS.time))+geom_point(aes(col=Status))+theme_bw()+
    scale_colour_manual(values =cols)+
    labs(x=xlab,y=b.labs)+
    geom_vline(xintercept=sum(fp_dat$Risktype=="Low"),colour="black", linetype="dotted",size=0.8)
  p2
  p.merge=mg_merge_plot(p1,p2,nrow=2,ncol=1,labels = labs)
  return(p.merge)
}
get_riskscore<-function(dat,os,os.time,step=T,direction=c("both", "backward", "forward")[1]){
  tcga_dat1 <- cbind(time=os.time,
                     status=os,
                     dat)
  tcga_dat1=crbind2DataFrame(tcga_dat1)
  colnames(tcga_dat1)=gsub('-','__',colnames(tcga_dat1))
  gene111=gsub('-','__',colnames(dat))
  fmla <- as.formula(paste0("Surv(time, status) ~"
                            ,paste0(gene111,collapse = '+')))
  cox <- coxph(fmla, data =as.data.frame(tcga_dat1))
  if(step==T){
    cox1 <- step(cox,direction =direction)
  }else{
    cox1=cox
  }
  lan <- coef(cox1)
  #round(lan, 3)
  genes <- names(cox1$coefficients)
  mult_results=paste0(round(lan, 3), '*', names(lan),collapse = '+')
  risk.tcga <- as.numeric(lan%*%as.matrix(t(tcga_dat1[,genes])))
  
  data_gene_score_final<-tcga_dat1
  data_gene_score_final$Samples<-rownames(data_gene_score_final)
  data_gene_score_final$riskscore=risk.tcga
  #data_gene_score_final$riskscorez=mosaic::zscore(risk.tcga)
  #最佳截断
  optimalCutoff <- survminer::surv_cutpoint(data.frame(time=data_gene_score_final$time/365,
                                                       event=data_gene_score_final$status,
                                                       risk=data_gene_score_final$riskscore), 
                                            time = "time", event = "event",variables = c("risk"))
  optimalCutoff=optimalCutoff$cutpoint$cutpoint[1]
  #print(optimalCutoff)
  #optimalCutoff=median(data_gene_score_final$riskscore)
  #optimalCutoff=0
  data_gene_score_final$Risk=ifelse(data_gene_score_final$riskscore>optimalCutoff,'High','Low')
  table(data_gene_score_final$Risk)
  data_gene_score_final$cutoff=optimalCutoff
  return(list(result=data_gene_score_final,module.gene=genes,model=mult_results))
}
alter_graphic <- function (graphic = c("rect", "point"), width = 1,
                           height = 1, horiz_margin = unit(1, "pt"), vertical_margin = unit(1,
                                                                                            "pt"), fill = "red", col = NA, pch = 16,
                           ...)
{
  graphic = match.arg(graphic)[1]
  if (graphic == "rect") {
    if (!is.numeric(width)) {
      stop_wrap("`width` should be nummeric.")
    }
    if (!is.numeric(height)) {
      stop_wrap("`height` should be nummeric.")
    }
    if (width != 1) {
      if (missing(horiz_margin)) {
        horiz_margin = unit(0, "pt")
      }
    }
    if (height != 1) {
      if (missing(vertical_margin)) {
        vertical_margin = unit(0, "pt")
      }
    }
    fun = function(x, y, w, h) {
      w = w * width
      h = h * height
      grid.rect(x, y, w - horiz_margin * 2, h - vertical_margin *
                  2, gp = gpar(fill = fill, col = col, ...))
    }
  }
  else if (graphic == "point") {
    fun = function(x, y, w, h) {
      grid.points(x, y, pch = pch, gp = gpar(fill = fill,
                                             col = col, ...))
    }
  }
  return(fun)
}
mg_violin_1=function(data,xangle=0,ylab='value',xlab='',leg.title='Group',test_method='anova',cmp_test_method='t.test',
                     legend.pos='r',melt=F,jitter=T,ylim=NULL,show_compare=NULL,point_size=NULL,group_col){
  library(ggplot2)
  if(is.null(ylim)){
    
  }
  if(melt){
    data_m=data
    colnames(data_m)=c('Group','value')
  }else{
    data_m=reshape2::melt(data)
    colnames(data_m)=c('Group','value')
  }
  if(!is.null(ylim)){
    data_m$value[data_m$value>ylim[2]]<-NA
  }
  data_m=data_m[which(!is.na(data_m[,1])),]
  if(xangle==0){
    tx=element_text(colour="black",family="Times")
  }else{
    tx=element_text(angle=xangle,hjust = 1,colour="black",family="Times")
  }
  
  pos='right'
  if(is.null(legend.pos)){
    pos='none'
  }else if(legend.pos=='tr'){
    pos=c(1,1)
  }else if(legend.pos=='br'){
    pos=c(1,0)
  }else if(legend.pos=='tl'){
    pos=c(0,1)
  }else if(legend.pos=='bl'){
    pos=c(0,0)
  }else if(legend.pos=='t'){
    pos='top'
  }else if(legend.pos=='r'){
    pos='right'
  }else if(legend.pos=='b'){
    pos='bottom'
  }
  uni.group=unique(data_m[,1])
  ct=length(uni.group)
  
  p1<-ggplot(data_m,aes(x=Group,y=value))+geom_violin(alpha=0.7)
  # if(ct<=4){
  #   p1=p1+ggsci::scale_fill_lancet()
  # }else if(ct<=10){
  #   p1=p1+ggsci::scale_fill_npg(name=leg.title)
  # }else if(ct<=20){
  #   p1=p1+ggsci::scale_fill_d3(palette = "category20",name=leg.title)
  # }else if(ct<=30){
  #   cbPalette=c(ggsci::pal_npg("nrc", alpha = 0.6)(10),ggsci::pal_d3("category20", alpha = 0.6)(20))
  #   p1=p1+scale_fill_manual(values=cbPalette[1:ct])
  # }else if(ct<=38){
  #   cbPalette=c(ggsci::pal_lancet()(10)
  #               ,ggsci::pal_npg("nrc", alpha = 0.6)(10)
  #               ,ggsci::pal_d3("category20", alpha = 0.6)(20)
  #               ,ggsci::pal_nejm("default", alpha = 0.6)(8))
  #   p1=p1+scale_fill_manual(values=cbPalette[1:ct])
  # }
  p1=p1+scale_fill_manual(values=group_col)
  # if(jitter){
  #   if(is.null(point_size)){
  #     p1<-p1+geom_jitter(alpha=0.3,col='black',show.legend=FALSE,width = 0.2)
  #   }else{
  #     p1<-p1+geom_jitter(alpha=0.3,col='black',show.legend=FALSE,width = 0.2,size=point_size)
  #   }
  # }
  # 
  p1=p1+theme_bw()+geom_boxplot(width=0.2,aes(fill=Group),outlier.shape = NA)
  p1=p1+theme(axis.text.x=tx, 
              axis.text.y=element_text(family="Times",face="plain"), 
              axis.title.y=element_text(family="Times",face="plain"), 
              legend.text=element_text(face="plain", family="Times", colour="black" 
              ),
              legend.title=element_text(face="plain", family="Times", colour="black"
              ),
              legend.justification=pos, legend.position=pos
              ,legend.background = element_rect(fill = NA, colour = NA)
  )+ylab(ylab)+xlab(xlab)
  til=''
  if(test_method=='anova'){
    if(length(unique(data_m[,1]))<3){
      x1=data_m[,2][which(data_m[,1]==unique(data_m[,1])[1])]
      x2=data_m[,2][which(data_m[,1]==unique(data_m[,1])[2])]
      pv=t.test(x1,x2)$p.value 
      til=paste0('t-tests p=',signif(pv,2))
    }else{
      fit <- aov(value~Group, data = data_m)
      pv=summary(fit)[[1]][5][[1]]
      fv=summary(fit)[[1]][4][[1]]
      til=paste0('ANOVA tests p=',signif(pv,2))
    }
  }else{
    if(length(unique(data_m[,1]))<3){
      x1=data_m[,2][which(data_m[,1]==unique(data_m[,1])[1])]
      x2=data_m[,2][which(data_m[,1]==unique(data_m[,1])[2])]
      pv=wilcox.test(x1,x2)$p.value 
      til=paste0('wilcox.tests p=',signif(pv,2))
    }else{
      fit=kruskal.test(value~Group, data = data_m)
      pv=fit$p.value
      til=paste0('Kruskal-Wallis test p=',signif(pv,2))
    }
  }
  p1=p1+ggtitle(til) 
  if(!is.null(ylim)){
    p1=p1+ylim(ylim)
  }
  if(is.null(show_compare)){
    if(length(uni.group)>5){
      show_compare=F
    }else{
      show_compare=T
    }
  }
  if(show_compare){
    comps=list()
    for(i in 1:(length(uni.group)-1)){
      for(j in (i+1):length(uni.group)){
        comps=c(comps,list(c(uni.group[i],uni.group[j])))
      }
    }
    p1=p1+ggpubr::stat_compare_means(comparisons = comps,method = cmp_test_method,label= "p.signif")
  }
  return(p1)
}
alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"),
              gp = gpar(fill = "#CCCCCC", col = NA))
  },
  
  Mutant = function(x, y, w, h) {
    grid.rect(x, y, w-unit(3, "pt"), h-unit(2, "pt"),
              gp = gpar(fill = col["Mutant"], col = NA))
  }
)

alter_fun = list(
  background = alter_graphic("rect", fill = "#CCCCCC"),
  Mutant = alter_graphic("rect", fill = col["Mutant"])
)
heatmap_legend_param = list(title = "Alterations",
                            at = c("Mutant"),
                            labels = c("Mutant"))
my_riskplot=function(cli_dat,cols=c("red","blue"),xlab='Samples',
                     a.ylab="Risk score",b.labs="Survival time(year)",cutoff=0,labs=c('A','B')){
  #cli_dat=tcga.Risktype.cli
  cli.dat.order=cli_dat[order(cli_dat$Riskscore),c('OS.time','Status','Riskscore','Risktype')]
  fp_dat=data.frame(Samples=1:nrow(cli_dat),cli.dat.order)
  p1=ggplot(fp_dat,aes(x=Samples,y=Riskscore))+geom_point(aes(color=Risktype))+
    scale_colour_manual(values =cols)+
    theme_bw()+labs(x=xlab,y=a.ylab)+
    geom_hline(yintercept=cutoff,colour="black", linetype="dotted",size=0.8)+
    geom_vline(xintercept=sum(fp_dat$Risktype=="Low"),colour="black", linetype="dotted",size=0.8)
  p1
  p2=ggplot(fp_dat,aes(x=Samples,y=OS.time))+geom_point(aes(col=Status))+theme_bw()+
    scale_colour_manual(values =cols)+
    labs(x=xlab,y=b.labs)+
    geom_vline(xintercept=sum(fp_dat$Risktype=="Low"),colour="black", linetype="dotted",size=0.8)
  p2
  p.merge=mg_merge_plot(p1,p2,nrow=2,ncol=1,labels = labs)
  return(p.merge)
}
get_riskscore.lasso<-function(dat,os,os.time,labels=c('A','B')){
  library(glmnet)
  set.seed(2021)
  fit1=glmnet(as.matrix(dat)
              ,cbind(time=os.time,
                     status=os)
              ,family="cox"
              ,nlambda=100
              , alpha=1) 
  
  cv.fit<-cv.glmnet(as.matrix(dat)
                    ,cbind(time=os.time,
                           status=os)
                    ,family="cox"
                    ,nfolds = 10
                    ,nlambda=100
                    , alpha=1)
  sig.coef <- coefficients(cv.fit,s=cv.fit$lambda.min)[which(coefficients(cv.fit,s=cv.fit$lambda.min)[,1]!=0),1]
  #print(cv.fit$lambda.min)
  #length(names(sig.coef))
  #10
  mg_plot_lasso <- function(fit,cv_fit,lambda=NULL,show_text=T,figLabels=c('A','B')){
    if(is.null(lambda)){
      lmda=cv_fit$lambda.min
    }else{
      lmda=lambda
    }
    fit.coef=fit$beta[(apply(fit$beta,1,function(x){
      return(sum(x!=0))
    })>0),]
    
    fit.coef=as.matrix(fit.coef)
    colnames(fit.coef)=fit$lambda
    #fit$lambda==cv_fit$lambda
    library(ggplot2)
    dat=data.table::melt(t(as.matrix(fit.coef)))
    dat_z=dat[which(dat$value==0),]
    dat=dat[which(dat$value!=0),]
    dat.sv=rbind()
    for (u in unique(dat_z[,2])) {
      t.z=dat_z[which(dat_z[,2]==u),1]
      t.zx=max(t.z)
      dat.sv=rbind(dat.sv,c(t.zx,u,0))
      t.zn=min(t.z)
      if(t.zx!=t.zn){
        dat.sv=rbind(dat.sv,c(t.zn,u,0))
      }
    }
    colnames(dat.sv)=colnames(dat_z)
    #dat_z=dat_z[dat_z[,2]%in%names(which(fit.coef[,which(fit$lambda==lmda)]!=0)),]
    dat=crbind2DataFrame(rbind(dat,dat.sv))
    mn=min(-log(dat$Var1))
    mx=max(-log(dat$Var1))
    if(show_text){
      mx=(mx-mn)*0.1+mx
    }
    p=ggplot(dat, aes(x=-log(Var1), y=value,colour=Var2))+geom_line()+theme_bw()+theme(legend.position = "none")
    p=p+coord_cartesian(xlim=c(mn, mx))+xlab('-ln(lambda)')+ylab('Coefficients')
    if(show_text){
      fl=fit.coef[which(fit.coef[,which(fit$lambda==lmda)]!=0),ncol(fit.coef)]
      for_label=data.frame(Var1=rep(min(dat$Var1),length(fl)),Var2=names(fl),value=fl)
      p=p+ggrepel::geom_label_repel(
        aes(label = Var2,color=Var2),
        data = for_label,hjust = 0
      )
    }
    p=p+geom_vline(aes(xintercept=-log(lmda)), colour="#BB0000", linetype="dashed")
    p=p+annotate('text',x=-log(lmda),y=min(dat[,3]),label=paste0('lambda=',round(lmda,4)))
    tgc=data.frame(lambda=cv_fit$lambda,cvm=cv_fit$cvm,cvup=cv_fit$cvup,cvlo=cv_fit$cvlo,cvsd=cv_fit$cvsd
                   ,col=ifelse(cv_fit$lambda>=cv_fit$lambda.min&cv_fit$lambda<=cv_fit$lambda.1se,ifelse(cv_fit$lambda==lmda,'A','C'),'B'))
    p1=ggplot(tgc, aes(x=log(lambda), y=cvm)) + xlab('ln(lambda)')+ ylab('Partial Likelihood Deviance')+
      geom_errorbar(aes(ymin=cvm-cvsd, ymax=cvm+cvsd)) +
      geom_point(aes(colour=col))
    p1=p1+theme_bw()+theme(legend.position = "none")
    gal=ggpubr::ggarrange(p,p1, ncol = 2, nrow = 1
                          #,align = "hv"
                          ,labels = figLabels)
    return(gal)
  }
  lasso.pdf <- mg_plot_lasso(fit1,
                             cv.fit,
                             show_text=T,
                             figLabels=labels)
  return(list(lasso.gene=names(sig.coef),lambda.min=cv.fit$lambda.min,plot=lasso.pdf))
}
ggplotKMCox=function(dat,title='Groups',labs=NULL,add_text=NULL, col = mycolor){
  library(ggplot2)
  colnames(dat)=c('time','status','groups')
  #sdf<-survdiff(Surv(time,status) ~ groups,data=dat)
  #print((sdf))
  #summary(sdf)
  #p<-pchisq(sdf$chisq,length(sdf$n)-1,lower.tail=FALSE)
  sf<-survfit(Surv(time,status) ~ groups,data=dat)
  surv=survminer::ggsurvplot(sf, data = dat, palette = col, #jco palette 
                             pval = TRUE, surv.median.line='hv',
                             linetype = "strata"
                             ,conf.int = F
                             # ,conf.int.style ='step'
                             , pval.coord=c(0, 0.2), #Add p-value 
                             risk.table = TRUE, 
                             legend.title = title
                             ,legend.labs = labs
  )
  p1=surv$plot+theme_bw()+theme(axis.text.y=element_text(family="Times",face="plain")
                                ,axis.text.x=element_blank()
                                ,axis.title.x=element_blank()
                                ,plot.margin=unit(c(0.2, 0.2, 0, 0.1), "inches")
                                #,axis.title.y=element_blank()
                                ,legend.position=c(1,1), legend.justification=c(1,1)
                                ,legend.background = element_rect(fill = NA, colour = NA)
                                ,legend.title = element_text(family="Times",face="plain")
                                ,legend.text = element_text(family="Times",face="plain"))
  #p1=p1+text()
  #tms=data.frame(Group=tms.gp,value=tms.tps,Attribute=rep(data_m[1,1],length(tms.gp))
  #               ,ymax=rep(max(ylim),length(tms.gp)))
  #p4=p4+geom_text(data=tms,aes(x=Group, y=ymax, label=value),color="black")
  if(!is.null(add_text)){
    text.tb=surv$data.survplot[1,]
    text.tb[1,1]=0
    text.tb[1,5]=0
    text.tb$Text=add_text
    p1=p1+geom_text(data=text.tb,aes(x=time, y=surv, label=Text),color="black",hjust =0)
  }
  
  p2=surv$table+theme_bw()+theme(axis.text.y=element_text(family="Times",face="plain")
                                 #,axis.text.x=element_blank()
                                 #,axis.title.x=element_blank()
                                 #,axis.title.y=element_blank()
                                 ,plot.margin=unit(c(0, 0.2, 0.2, 0.1), "inches")
                                 ,plot.title=element_blank()
                                 ,legend.position=c(1,1), legend.justification=c(1,1)
                                 #,legend.background = element_rect(fill = NA, colour = NA)
                                 ,legend.title = element_text(family="Times",face="plain")
                                 ,legend.text = element_text(family="Times",face="plain"))
  
  g2=ggpubr::ggarrange(p1,p2, ncol = 1, nrow = 2,heights = c(1,0.3),align = "v")
  return(g2)
}
mg_violin <- function(data,xangle=0,ylab='value',xlab='',leg.title='Group',test_method='anova',cmp_test_method='t.test',legend.pos='r',melt=F,jitter=T,ylim=NULL,show_compare=NULL,point_size=NULL,mycolor = mycolor){
  library(ggplot2)
  if(is.null(ylim)){
    
  }
  if(melt){
    data_m=data
    colnames(data_m)=c('Group','value')
  }else{
    data_m=reshape2::melt(data)
    colnames(data_m)=c('Group','value')
  }
  if(!is.null(ylim)){
    data_m$value[data_m$value>ylim[2]]<-NA
  }
  data_m=data_m[which(!is.na(data_m[,1])),]
  if(xangle==0){
    tx=element_text(colour="black",family="Times")
  }else{
    tx=element_text(angle=xangle,hjust = 1,colour="black",family="Times")
  }
  
  pos='right'
  if(is.null(legend.pos)){
    pos='none'
  }else if(legend.pos=='tr'){
    pos=c(1,1)
  }else if(legend.pos=='br'){
    pos=c(1,0)
  }else if(legend.pos=='tl'){
    pos=c(0,1)
  }else if(legend.pos=='bl'){
    pos=c(0,0)
  }else if(legend.pos=='t'){
    pos='top'
  }else if(legend.pos=='r'){
    pos='right'
  }else if(legend.pos=='b'){
    pos='bottom'
  }
  uni.group=unique(data_m[,1])
  ct=length(uni.group)
  
  p1<-ggplot(data_m,aes(x=Group,y=value))+geom_violin(alpha=1)+scale_color_manual(values = mycolor)
  if(ct<=4){
    p1=p1+ggsci::scale_fill_d3()
  }else if(ct<=10){
    p1=p1+ggsci::scale_fill_d3(name=leg.title)
  }else if(ct<=20){
    p1=p1+ggsci::scale_fill_d3(palette = "category20",name=leg.title)
  }else if(ct<=30){
    cbPalette=c(ggsci::pal_d3("nrc", alpha = 0.6)(10),ggsci::pal_d3("category20", alpha = 0.6)(20))
    p1=p1+scale_fill_manual(values=cbPalette[1:ct])
  }else if(ct<=38){
    cbPalette=c(ggsci::pal_jco()(10)
                ,ggsci::pal_d3("nrc", alpha = 0.6)(10)
                ,ggsci::pal_d3("category20", alpha = 0.6)(20)
                ,ggsci::pal_d3("default", alpha = 0.6)(8))
    p1=p1+scale_fill_manual(values=cbPalette[1:ct])
  }
  
  if(jitter){
    if(is.null(point_size)){
      p1<-p1+geom_jitter(alpha=0.3,col='black',show.legend=FALSE,width = 0.2)
    }else{
      p1<-p1+geom_jitter(alpha=0.3,col='black',show.legend=FALSE,width = 0.2,size=point_size)
    }
  }
  
  p1=p1+theme_bw()+geom_boxplot(width=0.2,aes(fill=Group),outlier.shape = NA)
  p1=p1+theme(axis.text.x=tx, 
              axis.text.y=element_text(family="Times",face="plain"), 
              axis.title.y=element_text(family="Times",face="plain"), 
              legend.text=element_text(face="plain", family="Times", colour="black" 
              ),
              legend.title=element_text(face="plain", family="Times", colour="black"
              ),
              legend.justification=pos, legend.position=pos
              ,legend.background = element_rect(fill = NA, colour = NA)
  )+ylab(ylab)+xlab(xlab)
  til=''
  if(test_method=='anova'){
    if(length(unique(data_m[,1]))<3){
      x1=data_m[,2][which(data_m[,1]==unique(data_m[,1])[1])]
      x2=data_m[,2][which(data_m[,1]==unique(data_m[,1])[2])]
      pv=t.test(x1,x2)$p.value 
      til=paste0('t-tests p=',signif(pv,2))
    }else{
      fit <- aov(value~Group, data = data_m)
      pv=summary(fit)[[1]][5][[1]]
      fv=summary(fit)[[1]][4][[1]]
      til=paste0('ANOVA tests p=',signif(pv,2))
    }
  }else{
    if(length(unique(data_m[,1]))<3){
      x1=data_m[,2][which(data_m[,1]==unique(data_m[,1])[1])]
      x2=data_m[,2][which(data_m[,1]==unique(data_m[,1])[2])]
      pv=wilcox.test(x1,x2)$p.value 
      til=paste0('wilcox.tests p=',signif(pv,2))
    }else{
      fit=kruskal.test(value~Group, data = data_m)
      pv=fit$p.value
      til=paste0('Kruskal-Wallis test p=',signif(pv,2))
    }
  }
  p1=p1+ggtitle(til) 
  if(!is.null(ylim)){
    p1=p1+ylim(ylim)
  }
  if(is.null(show_compare)){
    if(length(uni.group)>5){
      show_compare=F
    }else{
      show_compare=T
    }
  }
  if(show_compare){
    comps=list()
    for(i in 1:(length(uni.group)-1)){
      for(j in (i+1):length(uni.group)){
        comps=c(comps,list(c(uni.group[i],uni.group[j])))
      }
    }
    p1=p1+ggpubr::stat_compare_means(comparisons = comps,method = cmp_test_method,label= "p.signif", step_increase = 0.0)
  }
  return(p1)
}
ggplotTimeROC=function(time,status,score,mks=c(1,3,5), col = mycolor){
  #time=g.os
  #status=g.ev
  #score=as.numeric(cpm.score)
  #cx=coxRun(data.frame(time,status,score))
  #if(cx[1]<=1){
  #  score=-1*score
  #}
  roc.tm=mg_surv_pROC(time,status,score,mks)
  print('roc.tm')
  print((roc.tm))
  library(survival)
  library(ggplot2)
  mks=mg_predict_time_ymd(time,mks)
  print(mks)  
  ROC.DSST=timeROC::timeROC(T=time,
                            delta=status
                            ,marker=score,
                            cause=1,weighting="marginal",
                            times=mks,
                            iid=TRUE)
  print(ROC.DSST)
  mks=mks[which(!is.na(ROC.DSST$AUC)&ROC.DSST$AUC>0)]
  print(mks)
  if(length(mks)>0){
    if(max(ROC.DSST$AUC)<0.5){
      score=-1*score
    }
    ROC.DSST=timeROC::timeROC(T=time,
                              delta=status
                              ,marker=score,
                              cause=1,weighting="marginal",
                              times=mks,
                              iid=TRUE)
    print(ROC.DSST$times)
    if(max(ROC.DSST$times)<20){
      lb=paste0(ROC.DSST$times,'-Years')
    }else if(max(ROC.DSST$times)<365){
      lb=paste0(round(ROC.DSST$times/12,0),'-Years')
    }else{
      lb=paste0(round(ROC.DSST$times/365,0),'-Years')
    }
    
    lbs=paste0(lb,',AUC=',round(ROC.DSST$AUC,2),',95%CI(',paste0(round(confint(ROC.DSST,level = 0.95,na.rm=T)$CI_AUC[,1]/100,2),'-',
                                                                 round(confint(ROC.DSST,level = 0.95,na.rm=T)$CI_AUC[,2]/100,2)),')')
    #roc.tm=ROC.DSST$times[which(ROC.DSST$times>0)]
    
    #p.dat=rbind()
    #for(i in which(ROC.DSST$times>0)){
    #los=lowess(ROC.DSST$FP[,i], y=ROC.DSST$TP[,i], f = 1/3, iter = 100)
    #los$x=c(0,los$x,1)
    #los$y=c(0,los$y,1)
    # p.dat=rbind(p.dat,data.frame(los$x, y=los$y,rep(lbs[i],length(los$y)),stringsAsFactors = F))
    #}
    
    p.dat=rbind()
    print(length(roc.tm))
    for(i in 1:length(roc.tm)){
      #print(i)
      r1=roc.tm[[i]]
      x1=1-r1$specificities
      y1=r1$sensitivities
      #print(cbind(1-r1$specificities,r1$sensitivities))
      nx1=unique(x1)
      ny1=c()
      for(x in unique(x1)){
        x.inds=which(x1==x)
        if(length(x.inds)>0&x<0.5){
          ny1=c(ny1,min(y1[x.inds]))
        }else if(length(x.inds)>0){
          ny1=c(ny1,max(y1[x.inds]))
        }else{
          ny1=c(ny1,y1[x.inds][1])
        }
      }
      #print(cbind(nx1,ny1))
      p.dat=rbind(p.dat,data.frame(x=nx1, y=ny1,rep(lbs[i],length(nx1)),stringsAsFactors = F))
    }
    colnames(p.dat)=c('V1','V2','Type')
    p.dat=as.data.frame(p.dat)
    
    p1=ggplot(p.dat, aes(x=V1,y=V2, fill=Type))
    p1=p1+geom_line(aes(colour=Type),lwd=1.1)+
      theme_bw()+
      xlab('False positive fraction')+
      ylab('True positive fraction') +
      scale_color_manual(values = col)
    #p1=p1+stat_smooth(aes(colour=Type),se = FALSE, size = 1)+theme_bw()+xlab('False positive fraction')+ylab('True positive fraction') 
    
    p1=p1+theme(axis.text.y=element_text(family="Times",face="plain"),axis.text.x=element_text(family="Times",face="plain")
                ,axis.title.x=element_text(family="Times",face="plain"),axis.title.y=element_text(family="Times",face="plain")
                ,plot.title=element_blank()
                ,plot.margin=unit(c(0.1, 0.1, 0.1, 0.1), "inches")
                ,legend.position=c(1,0)
                ,legend.justification=c(1,0)
                ,legend.background = element_rect(fill = NA, colour = NA)
                ,legend.title = element_text(family="Times",face="plain")
                ,legend.text = element_text(family="Times",face="plain"))
    return(p1)
  }else{
    return(mg_getplot_bank('No data plot by ROC!'))
  }
}
mg_PlotMutiBoxplot=function(data,group,group_cols='jco'
                            ,test_method=c('t.test','wilcox.test','paired_t.test','paired_wilcox.test','anova','kruskal.test')[1]
                            ,order=NULL,size=1,fill=F,outlier.shape=NA,yscale=c('none','log2','log10')[1]
                            ,xangle=45,ylab='Value',xlab='',box_span=0.7
                            ,orientation = c("vertical", "horizontal", "reverse")[1]
                            ,legend.pos=NULL,melt=F,ylim=NULL,binwidth=0.05
                            ,add=c("none", "dotplot", "jitter", "boxplot", "point", "mean"
                                   , "mean_se", "mean_sd", "mean_ci", "mean_range", "median"
                                   , "median_iqr", "median_mad", "median_range")[3]){
  paired=FALSE
  if(test_method=='paired_t.test'|test_method=='paired_wilcox.test'){
    test_method=gsub('paired_','',test_method)
    paired=TRUE
  }
  print(class(data))
  if(add=='jitter'){
    fill=F
  }
  
  library(ggplot2)
  if(!melt){
    #print(class(data))
    if(class(data)=='numeric'|class(data)=='integer'){
      data=as.numeric(data)
      vd1.sbs=data.frame(group,rep('Tag',length(group)),data)
      #print(vd1.sbs)
    }else{
      data=as.data.frame(data)
      data$ID=group
      vd1.sbs <- reshape2::melt(data, id.vars=c("ID"))
    }
    colnames(vd1.sbs)=c('category','type','Score')
    Data=vd1.sbs
  }else{
    vd1.sbs=data
    colnames(vd1.sbs)=c('category','type','Score')
    Data=vd1.sbs
  }
  
  #vd1.sbs[,2]=paste0('C',as.numeric(as.character(vd1.sbs[,2])))
  if(is.null(order)){
    order=unique(vd1.sbs[,2])
  }
  
  if(xangle==0){
    tx=element_text(colour="black",family="Times")
  }else{
    tx=element_text(angle=xangle,hjust = 1,colour="black",family="Times")
  }
  
  pos=c(0,0)
  if(is.null(legend.pos)){
    pos='none'
  }else if(legend.pos=='tr'){
    pos=c(1,1)
  }else if(legend.pos=='br'){
    pos=c(1,0)
  }else if(legend.pos=='tl'){
    pos=c(0,1)
  }else if(legend.pos=='bl'){
    pos=c(0,0)
  }else if(legend.pos=='top'){
    pos='top'
  }else if(legend.pos=='buttom'){
    pos='buttom'
  }else{
    pos='right'
  }
  print(pos)
  if(fill){
    p <- ggpubr::ggboxplot(vd1.sbs, x="type", y="Score", fill = "category", yscale = yscale
                           ,palette = group_cols,width = box_span,size = size,order = order,outlier.shape=outlier.shape
                           ,orientation=orientation,add=add,add.params = list(binwidth=binwidth)
                           ,short.panel.labs = T)#
  }else{
    p <- ggpubr::ggboxplot(vd1.sbs, x="type", y="Score", color = "category", yscale = yscale
                           ,palette = group_cols,width = box_span,size = size,order = order,outlier.shape=outlier.shape
                           ,orientation=orientation,add=add,add.params = list(binwidth=binwidth)
                           ,short.panel.labs = T)#
  }
  
  p=p+ggpubr::stat_compare_means(aes(group=category), label = "p.signif", method = test_method,paired=paired
                                 #,label.y = max(vd1.sbs[,1])
  )
  #p=p+ylab(ylab)+xlab(xlab)
  #p=p+theme(axis.text.x = element_text(angle = xangle, hjust = 1))
  p=p+theme_bw()+theme(axis.text.x=tx, #
                       axis.text.y=element_text(family="Times",face="plain"), #
                       axis.title.y=element_text(family="Times",face="plain"), #
                       #panel.border = element_blank(),axis.line = element_line(colour = "black"), #
                       legend.text=element_text(face="plain", family="Times", colour="black"  #
                       ),
                       legend.title=element_text(face="plain", family="Times", colour="black" #
                       ),
                       legend.justification=pos, legend.position=pos
                       ,legend.background = element_rect(fill = NA, colour = NA)
                       #,panel.grid.major = element_blank(),   #
                       #panel.grid.minor = element_blank()
  )+ylab(ylab)+xlab(xlab) #
  if(!is.null(ylim)){
    p=p+ylim(ylim)
  }
  return(p)
}

plotMutiBar=function(dat,ist=F,margin=T,xlb='',ylb='',lineCol='black',lineW=0.5,legTitle='Group',showValue=F,showLine=T,xangle=0,isAuto=T){
  library(ggplot2)
  #library(tidyverse)
  #library(reshape2)
  #library(optparse)
  if(ist){
    dat=t(dat)
  }
  lbc=colnames(dat)
  lbr=row.names(dat)
  bk_dat=dat
  if(margin){
    dat=dat%*%diag(1/c(apply(t(dat), 1, sum)))
  }
  row.names(dat)=paste0('R',1:(nrow(dat)))
  colnames(dat)=paste0('C',1:(ncol(dat)))
  row.names(bk_dat)=paste0('R',1:(nrow(bk_dat)))
  colnames(bk_dat)=paste0('C',1:(ncol(bk_dat)))
  #df=cbind(bg=paste0('R',1:nrow(dat)),dat)
  #colnames(df)=c('bg',paste0('C',1:(ncol(dat))))
  tp.dat=as.data.frame(cbind(bg=row.names(dat),dat))
  tp.dat[,1]=as.character(tp.dat[,1])
  for(i in 2:ncol(tp.dat)){
    tp.dat[,i]=as.numeric(as.character(tp.dat[,i]))
  }
  mt.df=reshape2::melt(tp.dat)
  colnames(mt.df)=c('bg','variable','value')
  
  pg=ggplot(mt.df, aes(x=variable, y=value, fill=bg))+geom_bar(stat = "identity", width=lineW, col=lineCol) + theme_bw()
  if(showLine){
    for (i in 2:(ncol(tp.dat)-1)) {
      tmp=tp.dat[order(tp.dat[,1],decreasing = T),]
      tmp[,i]=base::cumsum(tmp[,i])
      tmp[,i+1]=base::cumsum(tmp[,i+1])
      colnames(tmp)[c(i,i+1)]=c('STY','ED')
      tmp1=cbind(tmp,STX=rep(i-1+lineW/2,nrow(tmp))
                 ,EDX=rep(i-lineW/2,nrow(tmp)))
      pg=pg+geom_segment(data=tmp1,aes(x=STX, xend=EDX, y=STY, yend=ED))
    }
  }
  
  if(showValue){
    pg=pg+geom_text(data=mt.df,aes(label=sprintf("%0.2f", round(value, digits = 2))),position=position_stack(vjust=0.5))
  }
  pg=pg+scale_x_discrete(breaks = paste0('C',1:(ncol(dat))),label = lbc)
  pg=pg+labs(x=xlb, y=ylb)+ggsci::scale_fill_d3()+theme(legend.position = "bottom")
  pg=pg+ggsci::scale_fill_d3()+scale_fill_discrete(breaks = paste0('R',1:nrow(dat)),label = lbr,name=legTitle)
  if(xangle>0){
    pg=pg+theme(axis.text.x = element_text(angle = xangle, hjust = 1),legend.position = "bottom")
  }
  
  g.tb=matrix(0,nrow=ncol(dat),ncol=ncol(dat))
  for(i in 1:(ncol(dat))){
    for(j in 1:ncol(dat)){
      if(i!=j){
        g.tb[i,j]=round(-log10((chisq.test(bk_dat[,c(i,j)])$p.value)),2)
      }
    }
  }
  colnames(g.tb)=lbc
  row.names(g.tb)=lbc
  g.tb=reshape2::melt(g.tb) 
  colnames(g.tb)=c('A1','A2','A3')
  g.tb$A4=paste0(g.tb[,3],ifelse(g.tb[,3]>-log10(0.05),'(*)',''))
  stable.p=ggplot(g.tb, aes(A1, A2)) + geom_tile(aes(fill = A3),colour = "white") + theme_bw()+xlab('')+ylab('')+ scale_fill_gradient(low = "white",high = "steelblue")+geom_text(aes(x=A1,y=A2,label=A4))+theme(legend.position="none",axis.title.x=element_blank(),axis.text.x=element_blank())
  stable.p=stable.p+ggtitle('-log10(anova p value)')
  if(isAuto){
    g1=ggpubr::ggarrange(stable.p,pg, ncol = 1, nrow = 2,heights = c(0.5,1),align = "hv")
    return(g1)
  }else{
    return(list(Bar=pg,Table=stable.p))
  }
}

mg_plot_lasso <- function(fit,cv_fit,lambda=NULL,show_text=T,figLabels=c('A','B')){
  if(is.null(lambda)){
    lmda=cv_fit$lambda.min
  }else{
    lmda=lambda
  }
  fit.coef=fit$beta[(apply(fit$beta,1,function(x){
    return(sum(x!=0))
  })>0),]
  
  fit.coef=as.matrix(fit.coef)
  colnames(fit.coef)=fit$lambda
  #fit$lambda==cv_fit$lambda
  library(ggplot2)
  dat=data.table::melt(t(as.matrix(fit.coef)))
  dat_z=dat[which(dat$value==0),]
  dat=dat[which(dat$value!=0),]
  dat.sv=rbind()
  for (u in unique(dat_z[,2])) {
    t.z=dat_z[which(dat_z[,2]==u),1]
    t.zx=max(t.z)
    dat.sv=rbind(dat.sv,c(t.zx,u,0))
    t.zn=min(t.z)
    if(t.zx!=t.zn){
      dat.sv=rbind(dat.sv,c(t.zn,u,0))
    }
  }
  colnames(dat.sv)=colnames(dat_z)
  #dat_z=dat_z[dat_z[,2]%in%names(which(fit.coef[,which(fit$lambda==lmda)]!=0)),]
  dat=crbind2DataFrame(rbind(dat,dat.sv))
  mn=min(-log(dat$Var1))
  mx=max(-log(dat$Var1))
  if(show_text){
    mx=(mx-mn)*0.1+mx
  }
  p=ggplot(dat, aes(x=-log(Var1), y=value,colour=Var2))+geom_line()+theme_bw()+theme(legend.position = "none")
  p=p+coord_cartesian(xlim=c(mn, mx))+xlab('-ln(lambda)')+ylab('Coefficients')
  if(show_text){
    fl=fit.coef[which(fit.coef[,which(fit$lambda==lmda)]!=0),ncol(fit.coef)]
    for_label=data.frame(Var1=rep(min(dat$Var1),length(fl)),Var2=names(fl),value=fl)
    p=p+ggrepel::geom_label_repel(
      aes(label = Var2,color=Var2),
      data = for_label,hjust = 0
    )
  }
  p=p+geom_vline(aes(xintercept=-log(lmda)), colour="#BB0000", linetype="dashed")
  p=p+annotate('text',x=-log(lmda),y=min(dat[,3]),label=paste0('lambda=',round(lmda,4)))
  tgc=data.frame(lambda=cv_fit$lambda,cvm=cv_fit$cvm,cvup=cv_fit$cvup,cvlo=cv_fit$cvlo,cvsd=cv_fit$cvsd
                 ,col=ifelse(cv_fit$lambda>=cv_fit$lambda.min&cv_fit$lambda<=cv_fit$lambda.1se,ifelse(cv_fit$lambda==lmda,'A','C'),'B'))
  p1=ggplot(tgc, aes(x=log(lambda), y=cvm)) + xlab('ln(lambda)')+ ylab('Parial Likelihood Deviance')+
    geom_errorbar(aes(ymin=cvm-cvsd, ymax=cvm+cvsd)) +
    geom_point(aes(colour=col))
  p1=p1+theme_bw()+theme(legend.position = "none")
  gal=ggpubr::ggarrange(p,p1, ncol = 2, nrow = 1
                        #,align = "hv"
                        ,labels = figLabels)
  return(gal)
}

cor_point <- function(x,
                      y,
                      method='Pearson',
                      top_col='red',
                      right_col='blue',
                      ylab='y expression',
                      xlab='x expression'){
  #x=rnorm(200)
  #y=rnorm(200)
  library(ggplot2)
  data=data.frame(x,y)  
  colnames(data)=c('wt','mpg')
  til=''
  if(method=='Pearson'){
    cr=cor.test(x,y)
    p=cr$p.value
    r=cr$estimate
    til=paste0('Pearson\'s correlation\nR=',round(r,3),'\nP=',signif(p,3))
  }else{
    cr=cor.test(x,y,method = "spearman")
    p=cr$p.value
    r=cr$estimate
    til=paste0('spearman correlation\nR=',round(r,3),'\nP=',signif(p,3))
  }
  
  empty <- ggplot()+geom_point(aes(1,1), colour="white") +
    theme(                              
      plot.background = element_blank(), 
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(), 
      panel.border = element_blank(), 
      panel.background = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank()
      ,plot.margin=unit(c(0.1, 0.1, 0, 0), "inches")
    )
  empty=empty+geom_text(aes(x=1, y=1, label=til),color="black")
  
  plot_top <- ggplot(data, aes(wt, fill=top_col)) + 
    geom_density(alpha=.5,fill=top_col) +theme_bw()+ 
    theme(legend.position = "none",                           
          #plot.background = element_blank(), 
          #panel.grid.major = element_blank(), 
          #panel.grid.minor = element_blank(), 
          panel.border = element_blank(), 
          #panel.background = element_blank(),
          axis.title.x = element_blank(),
          #axis.title.y = element_blank(),
          axis.text.x = element_blank(),
          #axis.text.y = element_blank(),
          axis.ticks.x = element_blank()
          ,plot.margin=unit(c(0.1, 0, 0, 0.1), "inches")
    )
  
  plot_right <- ggplot(data, aes(mpg, fill=right_col)) + 
    geom_density(alpha=.5,fill=right_col) +coord_flip()+theme_bw()+ 
    theme(legend.position = "none",                           
          #plot.background = element_blank(), 
          #panel.grid.major = element_blank(), 
          #panel.grid.minor = element_blank(), 
          panel.border = element_blank(), 
          #panel.background = element_blank(),
          #axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          #axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()
          ,plot.margin=unit(c(0.01, 0.1, 0.1, 0), "inches")
    )
  #scale_fill_manual(values = c("orange", "purple")) + 
  
  p1=ggplot(data=data, aes(x=wt, y=mpg))+geom_point()+stat_smooth(method="lm")
  p1=p1+theme_bw()
  p1=p1+theme(axis.text.x=element_text(family="Times",face="plain"), #
              axis.text.y=element_text(family="Times",face="plain"), #
              axis.title.y=element_text(family="Times",face="plain"), #
              #panel.border = element_blank(),
              axis.line = element_line(colour = "black"), #
              legend.text=element_text(face="plain", family="Times", colour="black"),  #
              legend.title=element_text(face="plain", family="Times", colour="black"), #
              plot.margin=unit(c(0.01, 0.01, 0.1, 0.1), "inches")
              #,panel.grid.major = element_blank(),   #
              #panel.grid.minor = element_blank()
  )+ylab(ylab)+xlab(xlab)
  
  pg1=ggpubr::ggarrange(plot_top,p1, ncol = 1, nrow = 2,heights = c(0.3,1),align = "v")
  pg2=ggpubr::ggarrange(empty,plot_right, ncol = 1, nrow = 2,heights = c(0.3,1),align = "v")
  
  pgal=ggpubr::ggarrange(pg1,pg2, ncol = 2, nrow = 1,widths = c(1,0.3),align = "h")
  return(pgal)
}

plot_GSEA_By_nodes_wb <- function(parseGSEAResult,
                                  indexs=c(1,2,3),
                                  color = mycolor,
                                  TermNames=NULL,
                                  left=NULL,
                                  right=NULL){
  #parseGSEAResult=gsea.result.kegg.result
  library(ggplot2)
  if(is.null(parseGSEAResult$TEMPLATE)){
    if(is.null(left)){
      left='RankTop'
    }
    if(is.null(right)){
      right='RankBottom'
    }
  }
  if(is.null(left)){
    left=parseGSEAResult$TEMPLATE[1]
  }
  if(is.null(right)){
    right=parseGSEAResult$TEMPLATE[2]
  }
  if(!is.null(TermNames)){
    inds=match(TermNames,parseGSEAResult$EnrichTable[,1])
    inds=inds[!is.na(inds)]
    if(length(inds)==0){
      print(paste0(TermNames,' Not Found!'))
      return(NA)
    }
  }else{
    inds=indexs
    if(max(inds)>nrow(parseGSEAResult$EnrichTable)){
      print(paste0(inds,' out range!'))
      return(NA)
    }
  }
  #parseGSEAResult=GSE17705_GSEA
  g.rnk=parseGSEAResult$Rank
  
  all.info=rbind()
  all.dat=rbind()
  for(i in inds){
    node=parseGSEAResult$Nodes[[i]]
    es.values=c(0,as.numeric(unlist(strsplit(XML::xmlGetAttr(node,'ES_PROFILE'),' '))),0)
    hit.index=c(0,as.numeric(unlist(strsplit(XML::xmlGetAttr(node,'HIT_INDICES'),' '))),nrow(g.rnk))
    m.inds=which.max(abs(es.values))
    es=as.numeric(XML::xmlGetAttr(node,'ES'))
    np=as.numeric(XML::xmlGetAttr(node,'NP'))
    FDR=as.numeric(XML::xmlGetAttr(node,'FDR'))
    nes=as.numeric(XML::xmlGetAttr(node,'NES'))
    title=gsub('^gene_sets.gmt#','',XML::xmlGetAttr(node,'GENESET'))
    length(hit.index)
    all.dat=rbind(all.dat,data.frame(Index=hit.index,ES=es.values,Term=rep(title,length(es.values))
                                     ,Type=c('A',rep('V',length(es.values)-2),'A')))
    all.info=rbind(all.info,c(title,es.values[m.inds[1]],
                              hit.index[m.inds[1]],es,nes,np,FDR))
  }
  all.info=crbind2DataFrame(all.info)
  #all.info
  
  cbPalette=color
  
  all.dat$Colour=cbPalette[as.numeric(as.factor(all.dat$Term))]
  
  col_mp=unique(cbind(as.character(all.dat$Term),all.dat$Colour))[,2]
  names(col_mp)=unique(cbind(as.character(all.dat$Term),all.dat$Colour))[,1]
  glb=unique(unlist(lapply(strsplit(all.info[,1],'_'),function(x){return(x[1])})))
  if(length(glb)==1){
    #g.labels=paste0(gsub(paste0('^',glb,'_'),'',g.labels))
    desc=gsub(paste0('^',glb,'_'),'',all.info[,1])
  }else{
    desc=all.info[,1]
  }
  ndesc=c()
  for(de in desc){
    #de=desc[1]
    if(nchar(de)>50){
      d2=paste0(substr(de,0,47),'...')
      ndesc=c(ndesc,d2)
    }else{
      ndesc=c(ndesc,de)
    }
  }
  g.labels=paste0(ndesc,'\nES=',signif(all.info[,4],2),',NES=',signif(all.info[,5],2),',P=',signif(all.info[,6],2),',FDR=',signif(all.info[,7],2))[match(names(col_mp),all.info[,1])]
  
  #dat=data.frame(Index=hit.index,ES=es.values)
  p=ggplot(data=all.dat, aes(x=Index, y=ES)) +geom_line(aes(colour=Term))+xlim(0,nrow(g.rnk))
  if(length(glb)==1){
    p=p+labs(title=paste0('Enrichment plot ',glb,' terms'))+theme(plot.title = element_text(hjust = 0.5))
  }
  p=p+scale_colour_manual(values=col_mp
                          ,breaks = names(col_mp)
                          ,labels = g.labels
  )
  p=p+ylab('Enrichment score')+theme_bw()
  #p+guides(color = FALSE)
  p=p+ geom_segment(aes(x = 0, xend = nrow(g.rnk), y = 0, yend = 0)
                    , color="grey"
                    ,linetype="dashed")
  
  p=p+theme(
    legend.position='none'
    ,axis.title.x=element_blank()
    ,axis.text.x = element_blank(),axis.ticks.x = element_blank()
    ,plot.margin=unit(c(0.2, 0.2, 0, 0.1), "inches"))
  #p+geom_label()
  es.min=min(all.dat$ES)
  ymin=es.min-(max(all.dat$ES)-es.min)*0.1
  
  lh=(es.min-ymin)*0.7
  
  dt1=all.dat
  dt2=all.dat
  dt1$Height=rep(ymin,nrow(all.dat))-lh*(as.numeric(as.factor(all.dat$Term))-1)
  dt2$Height=rep(ymin+(es.min-ymin)*0.7,nrow(all.dat))-lh*(as.numeric(as.factor(all.dat$Term))-1)
  dt1=dt1[which(dt1$Type!='A'),]
  dt2=dt2[which(dt2$Type!='A'),]
  #dt=rbind(dt1)
  p1=ggplot()
  #es.text=rbind()
  for(g in unique(dt1$Term)){
    cl=unique(dt1$Colour[which(dt1$Term==g)])[1]
    p1=p1+geom_line(data = rbind(dt1[which(dt1$Term==g),],dt2[which(dt1$Term==g),])
                    ,aes(x=Index,y=Height,group=Index),col=cl)
    #h1=dt1$Height[which(dt1$Term==g)][1]
    #es.text=rbind(es.text,c(g,0,h1))
  }
  #es.text=crbind2DataFrame(es.text)
  #all.info$SX=es.text[match(all.info[,1],es.text[,1]),2]
  #all.info$SY=es.text[match(all.info[,1],es.text[,1]),3]
  
  #p1=p1+geom_text(data = all.info,aes(x=SX,y=SY,label = paste0('ES=',signif(V4,2),',NES=',signif(V5,2),',P=',signif(V6,2),',FDR=',signif(V7,2)))
  #              ,vjust =-0, hjust = 0)
  
  p1=p1+theme_bw()+theme(legend.position='none',axis.title.y=element_blank(),axis.text.y = element_blank()
                         ,axis.title.x=element_blank(),axis.text.x = element_blank()
                         ,axis.line.x = element_blank(),axis.ticks.x = element_blank(),axis.ticks.y = element_blank()
                         ,axis.line.x.bottom = element_blank()
                         ,plot.margin=unit(c(0, 0.2, 0, 0.1), "inches")
                         ,axis.line = element_blank()
  )
  
  #ggpubr::ggarrange(p,p1, ncol = 1, nrow = 2,heights = c(1,0.1*length(inds)),align = "v")
  
  p2=ggplot(data=data.frame(Risk=c(0,g.rnk$V2,0),Index=c(1,1:nrow(g.rnk),nrow(g.rnk))),aes(y=Risk,x=Index))+geom_line()+theme_bw()
  p2=p2+ geom_segment(aes(x = 0, xend = nrow(g.rnk), y = 0, yend = 0)
                      , color="grey"
                      ,linetype="dashed")
  p2=p2+theme(plot.margin=unit(c(0, 0.2, 0.2, 0.1), "inches"))+ylab('Rank')+xlab('Rank in ordered dataset')
  p2=p2+geom_text(data=data.frame(Xl=c(0),Yl=c(0)),aes(x=0,y=0,label = left),vjust =1, hjust = 0)+geom_text(data=data.frame(Xl=c(nrow(g.rnk)),Yl=c(0)),
                                                                                                            aes(x=nrow(g.rnk),y=0,label = right),vjust =0, hjust = 1)
  g.h=0.1*length(inds)
  if(g.h>0.8){
    g.h=0.8
  }
  gal=ggpubr::ggarrange(p,p1,p2, ncol = 1, nrow = 3,heights = c(1,g.h,0.6),align = "v",common.legend = TRUE,legend = "right")
  return(gal)
}

plot_GSEA_By_node_wb <- function(parseGSEAResult,
                                 index=1,
                                 col = mycolor,
                                 TermName=NULL,
                                 left=NULL,
                                 right=NULL){
  library(ggplot2)
  if(is.null(parseGSEAResult$TEMPLATE)){
    if(is.null(left)){
      left='RankTop'
    }
    if(is.null(right)){
      right='RankBottom'
    }
  }
  if(is.null(left)){
    left=parseGSEAResult$TEMPLATE[1]
  }
  if(is.null(right)){
    right=parseGSEAResult$TEMPLATE[2]
  }
  if(!is.null(TermName)){
    ind=which(parseGSEAResult$EnrichTable[,1]==TermName)
    if(length(ind)==0){
      print(paste0(TermName,' Not Found!'))
      return(NA)
    }else{
      ind=ind[1]
    }
  }else{
    ind=index
    if(ind>nrow(parseGSEAResult$EnrichTable)){
      print(paste0(ind,' out range!'))
      return(NA)
    }
  }
  node=parseGSEAResult$Nodes[[ind]]
  g.rnk=parseGSEAResult$Rank
  es.values=c(0,as.numeric(unlist(strsplit(XML::xmlGetAttr(node,'ES_PROFILE'),' '))),0)
  hit.index=c(0,as.numeric(unlist(strsplit(XML::xmlGetAttr(node,'HIT_INDICES'),' '))),nrow(g.rnk))
  es=as.numeric(XML::xmlGetAttr(node,'ES'))
  nes=as.numeric(XML::xmlGetAttr(node,'NES'))
  np=as.numeric(XML::xmlGetAttr(node,'NP'))
  FDR=as.numeric(XML::xmlGetAttr(node,'FDR'))
  
  dat=data.frame(Index=hit.index,ES=es.values)
  p=ggplot(data=dat, aes(x=Index, y=ES)) +geom_line(aes(colour='darkgreen',size=2))+xlim(0,nrow(g.rnk))+theme_bw()+ylab('Enrichment score')+labs(title=gsub('^gene_sets.gmt#','',XML::xmlGetAttr(node,'GENESET')))+theme(plot.title = element_text(hjust = 0.5))
  p=p+ geom_segment(aes(x = 0, xend = nrow(g.rnk), y = 0, yend = 0)
                    , color="grey"
                    ,linetype="dashed")
  
  #p+geom_text(data=data.frame(Xl=c(0),yl=c(min(es.values))),aes(x=0,y=min(es.values),label = paste0('ES=',signif(es,2),'\nNES=',signif(nes,2)
  #                                                              ,'\nP=',signif(np,2),'\nFDR=',signif(FDR,2)))
  #            ,vjust =0, hjust = 0)
  
  if(es<0){
    p=p+geom_text(data=data.frame(Xl=c(0),yl=c(min(es.values))),aes(x=0,y=min(es.values),label = paste0('ES=',signif(es,2),'\nNES=',signif(nes,2),'\nP=',signif(np,2),'\nFDR=',signif(FDR,2)))
                  ,vjust =0, hjust = 0)
  }else{
    p=p+geom_text(data=data.frame(Xl=c(0),yl=c(min(es.values))),aes(x=nrow(g.rnk),y=max(es.values),label = paste0('ES=',signif(es,2),'\nNES=',signif(nes,2),'\nP=',signif(np,2),'\nFDR=',signif(FDR,2)))
                  ,vjust =1, hjust = 1)
  }
  es.min=min(dat$ES)
  if(es.min>0){
    es.min=0
  }
  
  ymin=es.min-(max(dat$ES)-es.min)*0.1
  dt1=dat
  dt2=dat
  dt1$Height=rep(ymin,nrow(dat))
  dt2$Height=rep(ymin+(es.min-ymin)*0.7,nrow(dat))
  p1=p+geom_line(data = rbind(dt1,dt2),aes(x=Index,y=Height,group=Index))
  p1=p1+ggforce::geom_link(data=data.frame(x=c(0,nrow(g.rnk)),y=c(ymin,ymin),xend=c(nrow(g.rnk)/2,nrow(g.rnk)/2),yend=c(ymin,ymin))
                           ,aes(x=x,y=y,xend=xend,yend=yend
                                ,alpha=1-..index..
                                ,colour=col
                                ,size=50
                           ))
  p1=p1+theme(legend.position='none',axis.title.x=element_blank(),axis.text.x = element_blank())
  p1=p1+geom_text(data=data.frame(Xl=c(0),yl=c(ymin)),
                  aes(x=0,y=ymin,label = left),vjust =0.5, hjust = 0)+
    geom_text(data=data.frame(Xl=c(nrow(g.rnk)),yl=c(ymin)),
              aes(x=nrow(g.rnk),y=ymin,label = right),vjust =0.5, hjust = 1)
  p1=p1+theme(plot.margin=unit(c(0.2, 0.2, 0, 0.1), "inches"))
  
  p2=ggplot(data=data.frame(Risk=c(0,g.rnk$V2,0),Index=c(1,1:nrow(g.rnk),nrow(g.rnk))),
            aes(y=Risk,x=Index))+geom_line()+theme_bw()
  p2=p2+ geom_segment(aes(x = 0, xend = nrow(g.rnk), y = 0, yend = 0)
                      , color="grey"
                      ,linetype="dashed")
  p2=p2+theme(plot.margin=unit(c(0, 0.2, 0.2, 0.1), "inches"))+ylab('Rank')+xlab('Rank in ordered dataset')
  gal=ggpubr::ggarrange(p1,p2, ncol = 1, nrow = 2,heights = c(1,0.6),align = "v")
  return(gal)
}

calculate_mRNAsi <- function(exp_file, genes = c('')) {
  library(gelnet)
  library(dplyr)
  # 
  norm <- read.delim("d:/public/PCBC/syn2701943_rnaseq_norm_gene_symbol.tsv") %>% tibble::column_to_rownames( "tracking_id" ) %>% as.matrix
  sam_type <- read.delim("d:/public/PCBC/syn2701943_rnaseq_norm_samples_type.tsv", 
                         header = T,stringsAsFactors = F, 
                         check.names = F)
  same_genes1 <- intersect(rownames(exp_file), rownames(norm))
  if (length(genes) > 1) {
    same_genes <- intersect(same_genes1, genes)
  } else {
    same_genes <- same_genes1
  }
  print(setdiff(genes, same_genes1))
  norm <- norm[same_genes, ]
  mean_norm <- apply(norm, 1, mean)
  norm <- norm - mean_norm
  norm.tr <- norm[,sam_type[sam_type$type == 'SC',]$Samples]
  norm.bk <- norm[,sam_type[sam_type$type != 'SC',]$Samples]
  mm <- gelnet(t(norm.tr), NULL, 0, 1)
  
  exp_file <- exp_file[same_genes, ]
  mRNAsi <- apply(exp_file, 2, function(z) {cor(z, mm$w, 
                                                method="sp", 
                                                use="complete.obs")})
  mRNAsi <- mRNAsi - min(mRNAsi)
  mRNAsi <- mRNAsi / max(mRNAsi)
  mRNAsi <- as.data.frame(cbind(mRNAsi), stringsAsFactors = FALSE)
  mRNAsi$samples <- rownames(mRNAsi)
  return(list(mRNAsi_index = mRNAsi,train = mm$w))
}

wb_beeswarm_plot <- function(dat = NULL,
                             show_compare = T,
                             xlab = 'Groups',
                             ylab = '',
                             method = c('t.test', 'wilcox.test')[1],
                             col = mycolor,
                             leg.pos = c('top','left','right','bottom','none')[1],
                             title = NULL,
                             group = 'Cluster') {
  library(ggbeeswarm)
  colnames(dat) <- c('Cluster', 'Feature')
  
  
  p1 <- ggplot(dat, aes(Cluster, Feature, color = Cluster)) + geom_quasirandom(method = "frowney") +
    ggtitle(title) + scale_color_manual(values = col[1:length(unique(dat$Cluster))]) +
    xlab(xlab) + ylab(ylab) + guides(color=guide_legend(title = group)) + theme_bw() +
    theme(legend.position=leg.pos)
  
  
  if(show_compare){
    uni.group = as.character(unique(dat$Cluster))
    comps=list()
    for(i in 1:(length(uni.group)-1)){
      for(j in (i+1):length(uni.group)){
        comps=c(comps,list(c(uni.group[i],uni.group[j])))
      }
    }
    p1=p1+ggpubr::stat_compare_means(comparisons = comps,
                                     method = method,
                                     label= "p.signif", 
                                     step_increase = 0.0)
  }
  return(p1)
}

library(ComplexHeatmap)

library(scales)
library(ggsci)
mycolor <- pal_d3(alpha =1)(7)
show_col(mycolor)


mycolors.scrna <- sample(mg_colors, 30)
mycolors.scrna <- c(mycolor, mycolors.scrna)
show_col(mycolors.scrna)

library(data.table)
library(stringr)


############
#TCGA-LAML################
library(data.table)
library(stringr)

tcga_cli <- read.delim('00_origin_datas/PMC6066282-TCGA-CDR-clinical.txt',
                       header = T, check.names = F, stringsAsFactors = F)
tcga_cli <- tcga_cli[!duplicated(tcga_cli$bcr_patient_barcode), ]
colnames(tcga_cli)
tcga_cli <- subset(tcga_cli,
                   !is.na(OS.time) &
                     OS.time > 0 & 
                     !is.na(OS))
tcga_cli$bcr_patient_barcode[tcga_cli$type == 'LAML'] <- paste0(tcga_cli$bcr_patient_barcode[tcga_cli$type == 'LAML'], '-03')
tcga_cli$bcr_patient_barcode[tcga_cli$type != 'LAML'] <- paste0(tcga_cli$bcr_patient_barcode[tcga_cli$type != 'LAML'], '-01')
rownames(tcga_cli) <- tcga_cli$bcr_patient_barcode
colnames(tcga_cli)[1] <- 'Samples'
tcga_types <- as.character(unique(tcga_cli$type))


tcga_exp<-read.delim('00_origin_datas/LAML_Merge_RNA_seq_FPKM.txt',sep='\t',header = T,row.names = 1,check.names = F)
tcga_exp[1:4,1:4]
table(substr(colnames(tcga_exp),14,15))
tcga_exp <- exp_ensg2symbol(tcga_exp)


genecode=read.delim('00_origin_datas/GeneTag.genecode.v32.txt',sep='\t',header = T)
table(genecode$TYPE)
mrna_genecode=genecode[which(genecode$TYPE=='protein_coding'),]
range(tcga_exp)
dim(tcga_exp)
# 25019   151
tcga_exp=log2(tcga_exp[intersect(rownames(tcga_exp),mrna_genecode$SYMBOL),]+1)
dim(tcga_exp)
#18448   151
range(tcga_exp)
#0.0000 13.9604

com.samples <- intersect(colnames(tcga_exp),
                         tcga_cli$Samples)
tcga_exp <- tcga_exp[, com.samples]
dim(tcga_exp)
tcga_cli <- tcga_cli[com.samples, ]



############
library(stringr)


# GSE10358####
GSE10358_cli <- getGEOSampleData('GSE10358')
GSE10358_cli <- GSE10358_cli[, c("Acc", "os months (3.31.10)", "vital status")]
colnames(GSE10358_cli) <- c("Samples", "OS.time", "OS")
GSE10358_cli <- GSE10358_cli[GSE10358_cli$OS != 'NULL', ]
table(GSE10358_cli$OS)
GSE10358_cli$OS <- ifelse(GSE10358_cli$OS == 'Alive', 0, 1)
GSE10358_cli$OS.time <- as.numeric(GSE10358_cli$OS.time) * 30
rownames(GSE10358_cli) <- GSE10358_cli$Samples

GSE10358 <- getGEOExpData('GSE10358')

GSE10358_anno <- GSE10358$Anno$GPL570

GSE10358_exp <- GSE10358$Exp$GPL570_54675_Data_col1
GSE10358_exp <- log2(GSE10358_exp + 1)
boxplot(GSE10358_exp[, 1:5])

GSE10358_exp <- exp_probe2symbol_v2(GSE10358_exp,
                                    GSE10358_anno[, c(1, 11)])
GSE10358_exp <- GSE10358_exp[str_split_fixed(rownames(GSE10358_exp), ' /// ', 2)[, 2] == '', ]
boxplot(GSE10358_exp[, 1:5])
# rm(GSE10358_cli, GSE10358_anno, GSE10358_exp, GSE10358)
dim(GSE10358_exp )
#21655   304
# GSE106291#####
GSE106291_cli <- getGEOSampleData('GSE106291')
GSE106291_cli <- GSE106291_cli[, c("Acc", "overall survival (days)", "life status", "Title")]
colnames(GSE106291_cli) <- c("Samples", "OS.time", "OS", "Title")
table(GSE106291_cli$OS)
GSE106291_cli$OS <- ifelse(GSE106291_cli$OS == 'alive', 0, 1)
rownames(GSE106291_cli) <- GSE106291_cli$Samples

GSE106291_exp <- openxlsx::read.xlsx('00_origin_datas/GEO/GSE106291_Matrix_table.xlsx',
                                     sheet = 1)
GSE106291_exp <- na.omit(GSE106291_exp)
rownames(GSE106291_exp) <- GSE106291_exp$GEO_ID
GSE106291_exp <- GSE106291_exp[, -1]
boxplot(GSE106291_exp[, 1:5])
GSE106291_exp[1:5, 1:5]
intersect(colnames(GSE106291_exp), GSE106291_cli$Title)
GSE106291_exp <- GSE106291_exp[, GSE106291_cli$Title]

colnames(GSE106291_exp) <- GSE106291_cli$Samples
dim(GSE106291_exp )
#21402   250


#########
rownames_list <- list(

  GSE106291_exp = rownames(GSE106291_exp),
  GSE10358_exp = rownames(GSE10358_exp),
  tcga_exp = rownames(tcga_exp)
)


intersections_gene <- Reduce(intersect,rownames_list)
length(intersections_gene)
##############
dir.create('01_SLC_gene')
SLC.gc=read.csv("01_SLC_gene/SLC_gene.csv",header = T)
SLC.gc=SLC.gc[grep("Solute Carrier Family",SLC.gc$Description),]
write.csv(SLC.gc,file = "01_SLC_gene/SLC.gc_select.csv")
SLC.gc <- SLC.gc$Gene.Symbol
length(SLC.gc)
length(intersect(SLC.gc,rownames(tcga_exp)))
setdiff(SLC.gc,rownames(tcga_exp))

######
SLC=SLC.gc
SLC.score=t(ssGSEAScore_by_genes(gene.exp = tcga_exp,genes = SLC))
#SLC.score=data.frame(Samples=rownames(SLC.score),score=as.numeric(SLC.score[,1]))
colnames(SLC.score)='SLC.score'
writeMatrix(SLC.score,'01_SLC_gene/SLC.score.txt',row = T)

#02.WGCNA####
dir.create('02_WGCNA')
library(WGCNA)
allowWGCNAThreads(nThreads = 36)
enableWGCNAThreads(nThreads = 36)

my_mad <- function(x){mad(x,na.rm = TRUE)}
wgcna_exp=t(tcga_exp)
m.mad <- apply(wgcna_exp,2,my_mad)
dim(tcga_exp)
# 18448   130

tpm_T2 <- wgcna_exp[,which(m.mad >max( quantile(m.mad, probs=seq(0, 1, 0.25))[2],0.01))]
# dim(tpm_T2)
#  60 13836
#tpm_T2=tcga_exp[which(apply(tcga_exp,1,sd)>0.5),]
#tpm_T2=(2^tpm_T2-1)
range(tpm_T2)

pdf('02_WGCNA/1.pdf',width = 10,height = 10)
tpm_T2.power=mg_wgcna_get_power(tpm_T2)
dev.off()

# minModuleSize = 30,    
# mergeCutHeight = 0.25, 
tpm_T2.power$cutPower
tpm_T2.module=mg_WGCNA_getModule(tpm_T2,
                                 power = tpm_T2.power$cutPower,
                                 deepSplit=2,
                                 mergeCutHeight=0.3,
                                 minModuleSize=60)


table(tpm_T2.module$Modules[,2])
length(table(tpm_T2.module$Modules[,2]))
#10
write.csv(tpm_T2.module$Modules,file = "02_WGCNA/WGCNA_Modules.csv",row.names = T)
pdf('02_WGCNA/2.pdf',height = 5,width = 6)
plotDendroAndColors(tpm_T2.module$Tree, tpm_T2.module$Modules,
                    c("Dynamic Module",'Merged Module'),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()
#
MODULE.gene.num <- as.data.frame(table(tpm_T2.module$Modules[,2]))
write.table(MODULE.gene.num,file = "02_WGCNA/MODULE.gene.num.txt",sep = "\t",quote = F,row.names =F)
pdf('02_WGCNA/3.pdf',height = 6,width = 6)
mg_barplot_point(labels = names(table(tpm_T2.module$Modules[,2]))
                 ,values = as.numeric(table(tpm_T2.module$Modules[,2]))
                 ,point_sizes = 2
                 ,point_cols = names(table(tpm_T2.module$Modules[,2]))
                 ,xlab = 'Number of Genes',legend.pos = NULL)
dev.off()

#### 
# Calculate eigengenes
MEs = tpm_T2.module$MEs
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
pdf('02_WGCNA/4.pdf',height = 6,width = 12,onefile = T)
plot(METree, main = "Clustering of module eigengenes",xlab = "", sub = "")
dev.off()



SLC.score

tcga_cli_use <-cbind.data.frame(tcga_cli[tcga_cli$Samples,],SLC.score=SLC.score[tcga_cli$Samples,])
head(tcga_cli_use)
tcga_cli_use=tcga_cli_use[,-1]
colnames(tcga_cli_use)
tcga_cli_use.part=data.frame(SLC.score=tcga_cli_use[tcga_cli$Samples,c(31)])
str(tcga_cli_use.part)
tcga_cli_use.part=sapply(tcga_cli_use.part, function(x)as.numeric(as.factor(x)))


spms=tcga_cli_use.part
MEs_col<-tpm_T2.module$MEs
dim(MEs_col)
modTraitCor = cor(MEs_col[,rownames(MEDiss)[METree$order]]
                  , spms
                  ,use = 'pairwise.complete.obs')
modTraitP = corPvalueStudent(modTraitCor, dim(spms)[1])
textMatrix = paste(signif(modTraitCor, 2), " (", format(modTraitP,scientific =TRUE,digits = 3), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)
dim(textMatrix)

rownames(modTraitCor)=gsub("ME","",rownames(modTraitCor))
rownames(textMatrix)=gsub("ME","",rownames(textMatrix))
colnames(modTraitCor)

pdf('02_WGCNA/5.pdf',width = 4,height =6)
labeledHeatmap(Matrix = data.frame(modTraitCor),
               xLabels = colnames(modTraitCor),
               yLabels = rownames(modTraitCor),
               cex.lab = 1,
               ySymbols = colnames(t(modTraitCor)), colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = data.frame(textMatrix),
               setStdMargins = FALSE,
               cex.text = 0.8, zlim = c(-1,1),
               main = paste("Module-trait relationships"),
               xLabelsAngle = 0)
dev.off()


geneModuleMembership <- signedKME(tpm_T2
                                  , data.frame(tpm_T2.module$MEs)
                                  , outputColumnName = "")
head(geneModuleMembership)
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership)
                                           , nrow(tpm_T2.module$MEs)))

geneTraitSignificance <- as.data.frame(cor(tpm_T2
                                           , spms
                                           , use = 'pairwise.complete.obs'))

GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance)
                                           , nrow(spms)))

modNames<-colnames(geneModuleMembership)
modNames

module = "turquoise"
column = match(module, modNames)
column


moduleGenes<- (tpm_T2.module$Modules[,'mergedColors']==module)

tcga.wgcna.gene=c(names(which(moduleGenes)))
length(tcga.wgcna.gene)
# 1607
write.table(tcga.wgcna.gene,file = "02_WGCNA/tcga.wgcna.gene.csv",sep = "\t",quote = F,row.names = F)
pdf('02_WGCNA/S_model_turquoise.pdf',width = 6,height = 6)
verboseScatterplot(geneModuleMembership[moduleGenes, column],
                   geneTraitSignificance[moduleGenes, 1],
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
dev.off()
######
dir.create("03_gsea")


######

cox.gene.enrichment=mg_clusterProfiler(genes =tcga.wgcna.gene )
fig3a=list()
fig3a[[1]]=enrichplot::dotplot(cox.gene.enrichment$KEGG)+ggtitle('KEGG')
fig3a[[2]]=enrichplot::dotplot(cox.gene.enrichment$GO_BP)+ggtitle('Biological Process')
fig3a[[3]]=enrichplot::dotplot(cox.gene.enrichment$GO_CC)+ggtitle('Cellular Component')
fig3a[[4]]=enrichplot::dotplot(cox.gene.enrichment$GO_MF)+ggtitle('Molecular Function')

fig3=mg_merge_plot(mg_merge_plot(fig3a[[1]],fig3a[[2]],labels = c('A','B'),ncol=2),
                   mg_merge_plot(fig3a[[3]],fig3a[[4]],ncol=2,nrow=1,labels = LETTERS[3:4],heights = c(1,1),widths = c(1,2)),
                   nrow=2,heights = c(1,1))

savePDF('03_gsea/Fig3.pdf',fig3,height = 12,width = 18)


########
dir.create('04_model')
######
tcga.SLC.score.gene <- data.frame(SLC_SCORE=SLC.score[rownames(tcga_cli),1],t(tcga_exp[tcga.wgcna.gene,rownames(tcga_cli)]))
outTab <- data.frame()
SLC_score <- "SLC_SCORE"                                  
j=tcga.SLC.score.gene[,SLC_score]
y=j
for (gene in colnames(tcga.SLC.score.gene)) {
  #i=tcga.SLC.score.gene[,colnames(tcga.SLC.score.gene)[2]]
  i=tcga.SLC.score.gene[,gene]
  x=i
  corT=cor.test(x,y)
  
  z=lm(y~x)
  cor=corT$estimate
  cor=round(cor,4)
  pvalue=corT$p.value
  pval=signif(pvalue,4)
  pval=format(pval, scientific = TRUE)
  
  outTab=rbind(outTab,
               cbind(SLC_score=SLC_score,gene=gene,cor=cor,
                     pvalue=pval))
}

corMatrix <- outTab[order(outTab$cor,decreasing = F),]
if(is.factor(corMatrix$cor)) {
  # Convert factor to character, then to numeric
  corMatrix$cor <- as.numeric(as.character(corMatrix$cor))
} else {
  # If it's already a character or numeric, just convert to numeric
  corMatrix$cor <- as.numeric(corMatrix$cor)
}
corMatrix$pvalue <- as.numeric(as.character(corMatrix$pvalue ))
cor=0.3
cor.p=0.01
corMatrix <-corMatrix[which(corMatrix$cor>cor&corMatrix$pvalue<cor.p),]
dim(corMatrix)
corMatrix <- corMatrix[-1,]
write.table(corMatrix,file = "04_model/corTable-SLC_SCORE-gene.xls",quote = F,sep = "\t",row.names = F)

######
sig.gene.cox=cox_batch(dat = tcga_exp[intersect(corMatrix$gene,rownames(tcga_exp)),tcga_cli$Samples],
                       time = tcga_cli$OS.time,event = tcga_cli$OS)
sig.gene.cox
write.csv(sig.gene.cox,'04_model/sig.cox.csv')


table(sig.gene.cox$p.value<0.01)
# FALSE  TRUE 
# 566   161
pre.genes=rownames(sig.gene.cox[sig.gene.cox$p.value<0.01,])
length(pre.genes)#161
tcga_model_data <- cbind(tcga_cli[, c("OS.time", "OS")],
                         t(tcga_exp[pre.genes, tcga_cli$Samples]))
colnames(tcga_model_data) <- gsub('-', '_', colnames(tcga_model_data))

######
tcga.lasso=get_riskscore.lasso(dat = tcga_model_data[,-c(1:2)],
                               os = tcga_model_data$OS,
                               os.time = tcga_model_data$OS.time)
length(tcga.lasso$lasso.gene)#9
tcga.lasso$plot

######
fmla <- as.formula(paste0("Surv(OS.time, OS) ~"
                          ,paste0(tcga.lasso$lasso.gene,collapse = '+')))
cox <- coxph(fmla, data =as.data.frame(tcga_model_data))
cox=step(cox)

lan <- coef(cox)
lan
paste0(round(lan, 3), '*', names(lan),collapse = '+')
write.csv(data.frame(gene=names(lan),coef=as.numeric(lan)),'04_model/gene_coef.csv',row.names = F)
####
######
module.coxforest=ggforest(cox, data = tcga_model_data, 
                          main = "Hazardratio", fontsize =1.0, 
                          noDigits = 2)
module.coxforest
ggsave('04_model/gene_forest.pdf',module.coxforest,height = 4,width = 9)
#######
lan.dataframe <- as.data.frame(lan)

lan.dataframe$gene <- rownames(lan.dataframe) 
lan.dataframe$gene <- factor(lan.dataframe$gene,levels = rownames(lan.dataframe)[order(lan.dataframe$lan)])

lan.dataframe$color_group <- ifelse(lan.dataframe$lan > 0, "Positive", "Negative")
library(ggplot2)




p <- ggplot(lan.dataframe, aes(x=gene, y=lan,fill=color_group)) +
  geom_bar(stat="identity") +
  xlab("Gene Name") +
  ylab("Coefficient") +
  ggtitle("Gene Coefficients") +
  coord_flip() +
  scale_fill_manual(values = c("Positive" = "#96CEB4", "Negative" = "#CC9999")) +
  theme_bw()+
  guides(fill=FALSE)
p1 <- p+geom_text(aes(label=sprintf("%.3f", lan)), hjust=-0.2, size=3, color="black")
p1
ggsave('04_model/gene_ Coefficients.pdf',height = 4,width = 6)


#######
Risktype.col=c('#CCCC00',"#C8A1E0")

risk.tcga=as.numeric(lan%*%as.matrix(t(tcga_model_data[tcga_cli$Samples,names(lan)])))
tcga.Risktype.cli=data.frame(tcga_cli,Riskscore=risk.tcga)

tcga.data.point <- surv_cutpoint(tcga.Risktype.cli, time = "OS.time", event = "OS",
                                 variables = 'Riskscore')
tcga.cutoff <- as.numeric(summary(tcga.data.point)[1])
tcga.cutoff
tcga.Risktype.cli$Risktype=ifelse(tcga.Risktype.cli$Riskscore>tcga.cutoff,'High','Low')

# #######
# tcga.Risktype.cli$Risktype=ifelse(tcga.Risktype.cli$Riskscore>median(tcga.Risktype.cli$Riskscore),'High','Low')

tcga.roc=ggplotTimeROC(tcga.Risktype.cli$OS.time,
                       tcga.Risktype.cli$OS,
                       tcga.Risktype.cli$Riskscore,mks = c(1,2,3,4,5))
tcga.roc
tcga.km.OS=ggsurvplot(fit=survfit(Surv(OS.time/365, OS) ~ Risktype,
                                  data = tcga.Risktype.cli),
                      data=tcga.Risktype.cli,
                      conf.int = T,pval = T,risk.table = T, 
                      fun = "pct",size = 1,surv.median.line = 'hv',
                      title='TCGA-LAML',legend.title='Risktype',
                      legend.labs = c('High','Low'),
                      linetype = c("solid", "dashed","strata")[1],
                      palette = Risktype.col,
                      ylab='Overall Survival(OS)',
                      legend=c(0.85,0.8),#标签位置
                      ggtheme = theme_bw(base_size = 12))
tcga.km.OS=mg_merge_plot(tcga.km.OS$plot,tcga.km.OS$table,nrow=2,heights = c(3,1),align = 'v')

tcga.km.OS

tcga.Risktype.cli$Status=ifelse(tcga.Risktype.cli$OS==0,'Alive','Dead')
tcga.model.p=my_riskplot(cli_dat = tcga.Risktype.cli,cols =Risktype.col,xlab = 'sample',
                         a.ylab = 'Riskscore',b.labs = 'Time(days)',cutoff = median(tcga.Risktype.cli$Riskscore),labs = '')




######
model.gene.df=data.frame(tcga_cli[,c('OS','OS.time')],t(tcga_exp[names(lan),tcga_cli$Samples]))
head(model.gene.df)
module.gene.km=list()
for (i in 1:length(names(lan))) {
  model.gene.df1=model.gene.df
  model.gene.df1$group=ifelse(model.gene.df1[,names(lan)[i]]>median(model.gene.df1[,names(lan)[i]]),'High','Low')
  module.gene.km[[i]]=ggsurvplot(fit=survfit(Surv(OS.time/365, OS) ~ group,
                                             data = model.gene.df1),
                                 data=model.gene.df1,
                                 conf.int = F,pval = T,risk.table = T, 
                                 fun = "pct",size = 1,surv.median.line = 'hv',
                                 title='TCGA',legend.title=names(lan)[i],
                                 # legend.labs = c('High','Low'),
                                 linetype = c("solid", "dashed","strata")[1],
                                 palette = Risktype.col,
                                 ylab='Overall Survival(OS)',
                                 legend=c(0.8,0.8),#标签位置
                                 ggtheme = theme_bw(base_size = 12))
  module.gene.km[[i]]=module.gene.km[[i]]$plot
}
tcga.module.km=mg_merge_plot(module.gene.km,ncol=4,nrow=1)

#tcga_risk_plot <- my_riskplot(tcga.risktype.cli,labs=c('D',''),cols = risktype.col)
###########
my_mutibarplot=function(df,xlab='group',leg.title='',cols=c('#134B70','#EF5A6F'))#pal_d3()(10)[5:6])
{
  prop.pval=round(chisq.test(df)$p.value,2)#round(-log10(chisq.test(df)$p.value),2)
  if( prop.pval<0.001)
    prop.pval='<0.001'
  df.prop=prop.table(df,margin=2)
  df.prop=reshape2::melt(df.prop)
  colnames(df.prop)<-c("type","group","Percentage")
  df.prop$Percentage<-round(df.prop$Percentage,digits=2)
  p=ggplot(df.prop,aes(x=group,y=Percentage,fill=type))+
    geom_bar(position = "fill",stat="identity")+
    scale_fill_manual(values = cols)+
    xlab(xlab)+labs(fill = leg.title,title = 'Chi-Squared Test',subtitle  =  paste0('pvalue  ',prop.pval))+
    theme_bw()+theme(text=element_text(family = 'Times'),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank())
  p
  return(p)
}
tcga.barplot=my_mutibarplot(df=table(tcga.risktype.cli$Status,tcga.risktype.cli$Risktype),xlab = 'Risktype',leg.title = 'Status')
tcga.barplot
######

###GSE106291####

GSE106291_model_data <- cbind(GSE106291_cli[, c("OS.time", "OS")],
                              t(GSE106291_exp[, GSE106291_cli$Samples]))
colnames(GSE106291_model_data) <- gsub('-', '_', colnames(GSE106291_model_data))
GSE106291_model_data[1:5, 1:5]


risk.GSE106291=as.numeric(lan%*%as.matrix(t(GSE106291_model_data[GSE106291_cli$Samples,names(lan)])))

GSE106291.Risktype.cli=data.frame(GSE106291_cli,Riskscore=risk.GSE106291)
# GSE106291.Risktype.cli$Risktype=ifelse(GSE106291.Risktype.cli$Riskscore>median(risk.GSE106291),'High','Low')
GSE106291.data.point <- surv_cutpoint(GSE106291.Risktype.cli, time = "OS.time", event = "OS",
                                     variables = 'Riskscore')
GSE106291.cutoff <- as.numeric(summary(GSE106291.data.point)[1])
GSE106291.cutoff
GSE106291.Risktype.cli$Risktype=ifelse(GSE106291.Risktype.cli$Riskscore>GSE106291.cutoff,'High','Low')
GSE106291.roc=ggplotTimeROC(GSE106291.Risktype.cli$OS.time,
                            GSE106291.Risktype.cli$OS,
                            GSE106291.Risktype.cli$Riskscore,mks = c(1,2,3,4,5))
GSE106291.roc
GSE106291.km=ggsurvplot(fit=survfit(Surv(OS.time/365, OS) ~ Risktype,
                                    data = GSE106291.Risktype.cli),
                        data=GSE106291.Risktype.cli,
                        conf.int = T,pval = T,risk.table = T,
                        fun = "pct",size = 1,surv.median.line = 'hv',
                        title='GSE106291',legend.title='Risktype',
                        legend.labs = c('High','Low'),
                        linetype = c("solid", "dashed","strata")[1],
                        palette = Risktype.col,ylab='Overall Survival(OS)',
                        legend=c(0.85,0.25),
                        ggtheme = theme_bw(base_size = 12))
GSE106291.km=mg_merge_plot(GSE106291.km$plot,GSE106291.km$table,nrow=2,heights = c(3,1),align = 'v')
GSE106291.km
GSE106291.risktype.cli$Status=ifelse(GSE106291.risktype.cli$OS==0,'alive','Dead')
GSE106291.barplot=my_mutibarplot(df=table(GSE106291.risktype.cli$Status,GSE106291.risktype.cli$Risktype),xlab = 'Risktype',leg.title = 'Status')
###GSE10358 roc 0.68####
GSE10358_model_data <- cbind(GSE10358_cli[, c("OS.time", "OS")],
                             t(GSE10358_exp[, GSE10358_cli$Samples]))
colnames(GSE10358_model_data) <- gsub('-', '_', colnames(GSE10358_model_data))
GSE10358_model_data[1:5, 1:5]


risk.GSE10358=as.numeric(lan%*%as.matrix(t(GSE10358_model_data[GSE10358_cli$Samples,names(lan)])))

GSE10358.Risktype.cli=data.frame(GSE10358_cli,Riskscore=risk.GSE10358)
#GSE10358.Risktype.cli$Risktype=ifelse(GSE10358.Risktype.cli$Riskscore>median(risk.GSE10358),'High','Low')
GSE10358.data.point <- surv_cutpoint(GSE10358.Risktype.cli, time = "OS.time", event = "OS",
                                     variables = 'Riskscore')
GSE10358.cutoff <- as.numeric(summary(GSE10358.data.point)[1])
GSE10358.cutoff
GSE10358.Risktype.cli$Risktype=ifelse(GSE10358.Risktype.cli$Riskscore>GSE10358.cutoff,'High','Low')
GSE10358.roc=ggplotTimeROC(GSE10358.Risktype.cli$OS.time,
                           GSE10358.Risktype.cli$OS,
                           GSE10358.Risktype.cli$Riskscore,mks = c(1,2,3,4,5))
GSE10358.roc
GSE10358.km=ggsurvplot(fit=survfit(Surv(OS.time/365, OS) ~ Risktype,
                                   data = GSE10358.Risktype.cli),
                       data=GSE10358.Risktype.cli,
                       conf.int = T,pval = T,risk.table = T,
                       fun = "pct",size = 1,surv.median.line = 'hv',
                       title='GSE10358',legend.title='Risktype',
                       legend.labs = c('High','Low'),
                       linetype = c("solid", "dashed","strata")[1],
                       palette = Risktype.col,ylab='Overall Survival(OS)',
                       legend=c(0.85,0.25),
                       ggtheme = theme_bw(base_size = 12))
GSE10358.km=mg_merge_plot(GSE10358.km$plot,GSE10358.km$table,nrow=2,heights = c(3,1),align = 'v')
GSE10358.km
GSE10358.risktype.cli$Status=ifelse(GSE10358.risktype.cli$OS==0,'alive','Dead')
GSE10358.barplot=my_mutibarplot(df=table(GSE10358.risktype.cli$Status,GSE10358.risktype.cli$Risktype),xlab = 'Risktype',leg.title = 'Status')
###########
pdf('04_model/Fig4.pdf',fig4,height = 20,width = 20)
fig4=mg_merge_plot(mg_merge_plot(tcga.lasso$plot,p1,labels = c('','C'),widths = c(1.7,1)),
                   mg_merge_plot(tcga.roc,tcga.km.OS,tcga.barplot,ncol=3,labels =LETTERS[4:6]),
                   mg_merge_plot(GSE10358.roc,GSE10358.km,GSE10358.barplot,ncol=3,labels = LETTERS[7:9]),
                   mg_merge_plot(GSE106291.roc,GSE106291.km,GSE106291.barplot,ncol=3,labels = LETTERS[10:12]),
                   nrow=4,heights = c(0.9,1,1,1))
dev.off()
ggsave('04_model/Fig4.pdf',fig4,height = 20,width = 20)


###########
dir.create('05_immune')

##5.1estimate####
tcga_estimate_score <- immu_estimate(exp = tcga_exp)
head(tcga_estimate_score)
dim(tcga_estimate_score)
dim(tcga.risktype.cli)
tcga_estimate_score_cli <- cbind(tcga.risktype.cli,
                                 tcga_estimate_score[tcga.risktype.cli$Samples, ])

tcga.risktype.cli
tcga_StromalScore_cor_plot <-  cor_point(x = tcga_estimate_score_cli$StromalScore,
                                         y = tcga_estimate_score_cli$Riskscore,
                                         top_col = mycolor[1],
                                         right_col = mycolor[2],
                                         xlab = 'StromalScore',
                                         ylab = 'RiskScore')
tcga_StromalScore_cor_plot

colnames(tcga_estimate_score_cli)
tcga_ImmuneScore_cor_plot <-  cor_point(x = tcga_estimate_score_cli$ImmuneScore,
                                        y = tcga_estimate_score_cli$Riskscore,
                                        top_col = mycolor[1],
                                        right_col = mycolor[2],
                                        xlab = 'ImmuneScore',
                                        ylab = 'RiskScore')
tcga_ImmuneScore_cor_plot

tcga_ESTIMATEScore_cor_plot <-  cor_point(x = tcga_estimate_score_cli$ESTIMATEScore,
                                          y = tcga_estimate_score_cli$Riskscore,
                                          top_col = mycolor[1],
                                          right_col = mycolor[2],
                                          xlab = 'ESTIMATEScore',
                                          ylab = 'RiskScore')
tcga_ESTIMATEScore_cor_plot



p5A<- cowplot::plot_grid(tcga_StromalScore_cor_plot,
                         tcga_ImmuneScore_cor_plot,
                         tcga_ESTIMATEScore_cor_plot,
                         ncol = 3,labels = c("A",'',''))
ggsave('05_immune/p5A.pdf',p5A,height = 4,width = 10)


##5.2ssgsea####
save(tcga_exp,file = "05_immune/tcga_exp.Rdata")
tcga.immu.ssgsea=immu_ssgsea(tcga_exp,isTCGA = T)
save(tcga.immu.ssgsea,file ='05_immune/tcga.immu.ssgsea.RDS')
tcga.immu.ssgsea <- readRDS("05_immune/tcga.immu.ssgsea.RDS")
fig5b=my_mutiboxplot(dat = tcga.immu.ssgsea[tcga.risktype.cli$Samples,],
                     group = tcga.risktype.cli$Risktype,
                     group_cols = risktype.col,
                     ylab = 'score')+ggtitle('ssGSEA')
ggsave('05_immune/immu.ssgsea.pdf',fig5b,height = 6,width = 6)

##5.3 ICGs####
tcga.icgs=immu_ICGs(tcga_exp)

icg.dat.RS=cbind(tcga.risktype.cli$Riskscore
                 ,tcga_model_data[tcga.risktype.cli$Samples,names(lan)]
                 ,tcga.icgs[tcga.risktype.cli$Samples,c('CTLA4','PDCD1','PDCD1LG2','LGALS9','CD80','CD28','HAVCR2')])
#c('CTLA4','PDCD1','PDCD1LG2','LGALS9','CD80','CD28','HAVCR2')
colnames(icg.dat.RS)[1]='Riskcsore'

icg_cor_res <- Hmisc::rcorr(as.matrix(icg.dat.RS),type = 'spearman')
icg_cor_res$P[is.na(icg_cor_res$P)] <- 0
icg_cor_res.p=icg_cor_res$P
icg_cor_res.p[1:5,1:5]
icg_cor_res.p<-ifelse(icg_cor_res.p<0.0001,'****',
                      ifelse(icg_cor_res.p<0.001,'***', 
                             ifelse(icg_cor_res.p<0.01,'**',
                                    ifelse(icg_cor_res.p<0.05,'*',''))))

pdf('05_immune/p5C.pdf',height = 6,width = 7,onefile = F)
pheatmap(icg_cor_res$r[-c(1:6),c(names(lan),'Riskcsore')],
         color = circlize::colorRamp2(c(-1, 0, 1), c('#FFDE4D', 'white', '#FF4C4C')),
         main="Heatmap", 
         display_numbers = icg_cor_res.p[-c(1:6),c(names(lan),'Riskcsore')],
         cluster_cols = F, 
         cluster_rows = F,
         show_rownames = T, 
         show_colnames = T,
         fontsize_row = 12, 
         fontsize_col = 16)
dev.off()



###########
dir.create('06_tide_drug')
#####################
tcga_tide_dat <- t(scale(t(tcga_exp),scale = F))
dim(tcga_tide_dat)
write.table(tcga_tide_dat,file = '06_tide_drug/tcga_tide_dat.txt',quote = F, sep = '\t')
tcga_tide_res<-read.csv('06_tide_drug/LAML_TIDE.csv',row.names = 1,stringsAsFactors = F)
head(tcga_tide_res)
tcga_tide_res=cbind(tcga_tide_res[tcga.risktype.cli$Samples,],tcga.risktype.cli)

tide_sel=c('TIDE','IFNG','Exclusion','Dysfunction','MDSC')
tcga_tide_list <- list()
for (fea in c("TIDE","Dysfunction","Exclusion","MDSC","CAF","TAM.M2")) {
  print(fea)
  tmp_plot <- mg_violin_1(data.frame(tcga.risktype.cli$Risktype
                                     ,tcga_tide_res[tcga.risktype.cli$Samples, fea])
                          ,melt = T
                          ,ylab = fea
                          ,jitter=T
                          ,group_col = risktype.col#pal_jco()(9)[1:2]
                          ,test_method = 'wilcox.test'
                          ,cmp_test_method = 'wilcox.test'
                          ,legend.pos = NULL
                          ,show_compare = T)
  tcga_tide_list[[fea]] <- tmp_plot
}
fig6a <- cowplot::plot_grid(plotlist = tcga_tide_list,
                            ncol = 6)
fig6a




############################
library(pRRophetic)
library(ggplot2)

# set.seed(12345)
# predictedPtype_Cisplatin <- pRRopheticPredict(as.matrix(tcga_exp)
#                                               , "Cisplatin"
#                                               , selection=1
#                                               ,dataset = "cgp2016")
# predictedPtype_Cisplatin <- data.frame(predictedPtype_Cisplatin)
# 
# tcga_durg_ic50_res <- predictedPtype_Cisplatin
# 
# drugs <- c("Cisplatin","Erlotinib","Rapamycin","Sunitinib","PHA-665752","MG-132","Paclitaxel","Cyclopamine","AZ628","Sorafenib","VX-680","Imatinib","TAE684","Crizotinib","Saracatinib","S-Trityl-L-cysteine","Z-LLNle-CHO","Dasatinib","GNF-2","CGP-60474","CGP-082996","A-770041","WH-4-023","WZ-1-84","BI-2536","BMS-509744","CMK","Pyrimethamine","JW-7-52-1","A-443654","GW843682X","MS-275","Parthenolide","KIN001-135","TGX221","Bortezomib","XMD8-85","Roscovitine","Salubrinal","Lapatinib","Vinorelbine","NSC-87877","QS11","CP466722","Midostaurin","Shikonin","AKT inhibitor VIII","Embelin","Bexarotene","Bleomycin","Phenformin")
# length(drugs)
# for (drug in drugs) {
#   print(drug)
#   set.seed(12345)
#   tmpic50 <- pRRopheticPredict(as.matrix(tcga_exp)
#                                , drug,tissueType = "digestive_system"
#                                , selection=1
#                                , dataset = "cgp2016")
#   tmpic50 <- data.frame(tmpic50)
#   colnames(tmpic50) <- drug
#   tcga_durg_ic50_res <- cbind(tcga_durg_ic50_res, tmpic50)
# }
# tcga_durg_ic50_res <- tcga_durg_ic50_res[, -1]
# save(tcga_durg_ic50_res,file='06_tide_drug/tcga_durg_ic50_res.Rdata')
load('06_tide_drug/tcga_durg_ic50_res.Rdata')
head(tcga_durg_ic50_res)
#tcga_durg_ic50_res=tcga_durg_ic50_res[,-1]
###
library(ggcorrplot)
library(psych)
IC50_RS_cor <- corr.test(x =tcga.risktype.cli$Riskscore,
                         y = tcga_durg_ic50_res[tcga.risktype.cli$Samples,],
                         method = "spearman",adjust = "BH",ci = F)


IC50_RS_cor_res=data.frame(drugs=colnames(tcga_durg_ic50_res))
IC50_RS_cor_res$cor<-as.numeric(IC50_RS_cor$r)
IC50_RS_cor_res$p.adj<-as.numeric(IC50_RS_cor$p.adj)
head(IC50_RS_cor_res)
table(IC50_RS_cor_res$p.adj<0.01)
# FALSE  TRUE 
# 63    23
IC50_RS_cor_res <- IC50_RS_cor_res[IC50_RS_cor_res$p.adj<0.01,]
IC50_RS_cor_res <- IC50_RS_cor_res [order(abs(IC50_RS_cor_res$cor),decreasing = T),]
drug_name <- IC50_RS_cor_res$drugs[c(1,4,5,6,10,14,17,20)]
###
drug.df=data.frame(tcga_durg_ic50_res[tcga.risktype.cli$Samples,drug_name],Risktype=tcga.risktype.cli$Risktype)
drug.df=melt(drug.df)
head(drug.df)
pdf('06_tide_drug/Fig6B.pdf',height = 4,width = 12)
dodge_width <-1.5

fig6b <- ggplot(drug.df, aes(x=variable, y=value, fill=Risktype)) +
  geom_violin(position = position_dodge(dodge_width), trim = FALSE,show.legend = T) +
  geom_boxplot(width = 0.3, position = position_dodge(dodge_width), outlier.shape = NA,show.legend = F) +
  scale_fill_manual(values = risktype.col) +
  facet_wrap(~variable, scales = 'free', nrow = 1, ncol = 8) +
  ggpubr::stat_compare_means(aes(group=Risktype), label = 'p.signif', method = 'wilcox.test') +
  ylab('IC50') + xlab('') 
#theme(legend.position = 'top')
dev.off()
fig6 <- mg_merge_plot(fig6a,fig6b,nrow = 2,common.legend = T,labels = c("A","B"))
ggsave('06_tide_drug/fig6.pdf',fig6,height = 8,width = 15)

#####
dir.create('07_risktype.mut')

########
tcga.maf=getTCGAMAFByCode('LAML')
tcga.risktype.use=tcga.risktype.cli[,c('Samples','Risktype')]
table(tcga.risktype.use$Risktype)
colnames(tcga.risktype.use)[1]='Tumor_Sample_Barcode'
tcga.risktype.use$Tumor_Sample_Barcode=substr(tcga.risktype.use$Tumor_Sample_Barcode,1,12)
tcga.risktype.use.high=tcga.risktype.use[which(tcga.risktype.use$Risktype=='High'),]
tcga.risktype.use.low=tcga.risktype.use[which(tcga.risktype.use$Risktype=='Low'),]

write.table(tcga.risktype.use.high,file='07_risktype.mut/tcga.risktype.use.high.txt')
write.table(tcga.risktype.use.low,file='07_risktype.mut/tcga.risktype.use.low.txt')

tcga.maf.high=subsetMaf(tcga.maf,tsb=intersect(tcga.maf@data$Tumor_Sample_Barcode,tcga.risktype.use.high$Tumor_Sample_Barcode))
tcga.maf.high<-read.maf(tcga.maf.high@data,isTCGA=T,clinicalData = '07_risktype.mut/tcga.risktype.use.high.txt')
tcga.maf.high@clinical.data

tcga.maf.low=subsetMaf(tcga.maf,tsb=intersect(tcga.maf@data$Tumor_Sample_Barcode,tcga.risktype.use.low$Tumor_Sample_Barcode))
tcga.maf.low<-read.maf(tcga.maf.low@data,isTCGA=T,clinicalData = '07_risktype.mut/tcga.risktype.use.low.txt')
tcga.maf.low@clinical.data

#######
pdf('07_risktype.mut/Fig7A.pdf',height = 4,width = 7,onefile = F)
oncoplot(maf = tcga.maf.high,top = 10,sortByAnnotation = T)
dev.off()
pdf('07_risktype.mut/Fig7B.pdf',height = 4,width = 7,onefile = F)
oncoplot(maf = tcga.maf.low,top = 10,sortByAnnotation = T)
dev.off()

tcga.tmb=mg_getTCGATMBByCode('LAML')
tcga.tmb$Sample=paste0(tcga.tmb$Sample,'-03')
tcga.tmb.ri=merge(tcga.tmb,tcga.risktype.cli,by.x='Sample',by.y='Samples')
rownames(tcga.tmb.ri)=tcga.tmb.ri$Sample
head(tcga.tmb.ri)
tcga.tmb.cli=tcga.tmb.ri
library(survival)
library(survminer)
tmb.cutoff<-surv_cutpoint(tcga.tmb.cli,
                          time="OS.time",
                          event="OS",
                          variables=c("TMB"))
summary(tmb.cutoff)
tcga.tmb.cli$type <- ifelse(tcga.tmb.cli$TMB > tmb.cutoff$cutpoint$cutpoint, 'TMB-High', 'TMB-Low')
tcga.tmb.cli$ri_tmb=rep('none',nrow(tcga.tmb.cli))
tcga.tmb.cli$ri_tmb[which(tcga.tmb.cli$Risktype=='High' & tcga.tmb.cli$type=='TMB-High')]='H-Risk & H-TMB'
tcga.tmb.cli$ri_tmb[which(tcga.tmb.cli$Risktype=='Low' & tcga.tmb.cli$type=='TMB-Low')]='L-Risk & L-TMB'
tcga.tmb.cli$ri_tmb[which(tcga.tmb.cli$Risktype=='Low' & tcga.tmb.cli$type=='TMB-High')]='L-Risk & H-TMB'
tcga.tmb.cli$ri_tmb[which(tcga.tmb.cli$Risktype=='High' & tcga.tmb.cli$type=='TMB-Low')]='H-Risk & L-TMB'
tcga.tmb.cli$ri_tmb[which(tcga.tmb.cli$ri_tmb=='none')]=NA
table(tcga.tmb.cli$ri_tmb)
pdf("07_risktype.mut/Fig7C.pdf",height = 6,width = 6)
fig7c=ggsurvplot(fit=survfit(Surv(OS.time, OS) ~ type,
                             data = data.frame(OS.time = tcga.tmb.cli$OS.time/365
                                               , OS = tcga.tmb.cli$OS
                                               , type=tcga.tmb.cli$ri_tmb)),
                 data=data.frame(OS.time = tcga.tmb.cli$OS.time/365
                                 , OS = tcga.tmb.cli$OS
                                 , type=tcga.tmb.cli$ri_tmb),
                 conf.int = T,pval = T,fun = "pct",risk.table = T, size = 0.7,surv.median.line = 'hv',
                 title='TMB & Risktype',ggtheme=theme_classic(),
                 linetype = c("solid", "dashed","strata")[1],
                 palette = pal_lancet()(9)[4:8],
                 #legend = c('top', 'bottom', 'left', 'right', 'none')[1],
                 legend = c(0.8,0.75), 
                 legend.title = "")#,legend.labs =  c('TMB-High','TMB-Low')
fig7c

dev.off()

################
library(progeny)
tcga.pathway.activ=progeny(as.matrix(tcga_exp),scale = T)
dim(tcga.pathway.activ)
range(tcga.pathway.activ)
mg_PlotMutiBoxplot(tcga.pathway.activ[tcga.risktype.cli$Samples,]
                   , group = tcga.risktype.cli$Risktype
                   , legend.pos = 'top'
                   #, group_cols = cluster.color
                   , test_method =  c('kruskal.test','wilcox.test','anova')[2]
                   , add = 'boxplot'
                   , ylab = 'score')


pathway_cor_RS=cbind.data.frame(Riskscore=tcga.risktype.cli$Riskscore,
                                tcga.pathway.activ[tcga.risktype.cli$Samples,])
cor_res <- Hmisc::rcorr(as.matrix(pathway_cor_RS),type = 'spearman')
cor_res$P[is.na(cor_res$P)] <- 0

pdf('07_risktype.mut/Fig7D.pdf',height = 7,width = 7)
fig7d <- corrplot(as.matrix(cor_res$r),
                  p.mat = as.matrix(cor_res$P),
                  mar = c(0,0,1,1),diag = F,
                  col=COL2('PRGn'),#diverging_hcl(100,palette = 'Green-Orange')
                  tl.srt = 90,tl.cex = 1,tl.col = 'black',tl.offset = 0.5,
                  cl.pos = c("b","r","n")[1],cl.align.text = 'l',cl.length = 5,
                  cl.ratio = 0.1,cl.cex = 0.8,
                  addgrid.col = 'white',
                  method = c("circle", "square", "ellipse", "number", "shade", "color", "pie")[6],
                  insig = 'label_sig',
                  sig.level=c(0.001,0.01,0.05),
                  pch.cex=1,is.corr=T,xpd=T)
fig7d
dev.off()

##############
dir.create('08_risktype.treatment')
library("IMvigor210CoreBiologies")
data(cds)
pheno<-pData(cds)
head(pheno)
exper_tpm=mg_get_immu_pd1_treament_exp()
exper_id=exper_tpm$tpm
exper_id$symbol=rownames(exper_id)
rownames(exper_id)<-exper_id$symbol
exper_id$symbol<-NULL
range(exper_id)
exper_id_use<-log2(exper_id+1)
dim(exper_id_use)
# 31085   348
range(exper_id_use)
#rownames(exper_id_use)=gsub('-','__',rownames(exper_id_use))
exper_id_use[1:5,1:5]


IMvigor210_model_data=data.frame(OS=pheno$censOS,OS.time = pheno$os,
                                 t(exper_id_use[intersect(names(lan),rownames(exper_id_use)),rownames(pheno)]))
head(IMvigor210_model_data)


imv210.module.risk=get_riskscore(dat = IMvigor210_model_data[,],
                                 os = IMvigor210_model_data$OS,
                                 os.time = IMvigor210_model_data$OS.time,
                                 step = F,direction = c("both", "backward", "forward")[1])
length(imv210.module.risk$module.gene)
imv210.module.risk$model

imv210.module.risk$result$Risktype=ifelse(imv210.module.risk$result$riskscore>median(imv210.module.risk$result$riskscore),'High','Low')


fig8a=ggsurvplot(fit = survfit(Surv(time,status)~Risktype,
                               data =imv210.module.risk$result),
                 data=imv210.module.risk$result,fun = "pct",
                 risk.table = T, title='IMvigor210',pval = T,conf.int = T,
                 legend = c(0.8,0.75),legend.labs = c("High","Low"),palette = risktype.col )

fig8a=mg_merge_plot(fig8a$plot,fig8a$table,nrow=2,heights = c(3,1),align = 'v')

imv.risk<-cbind.data.frame(IMvigor210_model_data[rownames(pheno),c('OS','OS.time')],
                           Riskscore=imv210.module.risk$result[rownames(pheno),]$riskscore,
                           Risktype=imv210.module.risk$result[rownames(pheno),]$Risktype,
                           data.frame(binaryResponse=pheno$binaryResponse,
                                      Response=pheno$`Best Confirmed Overall Response`,
                                      IC=pheno$`IC Level`,
                                      TC=pheno$`TC Level`,
                                      IP=pheno$`Immune phenotype`,
                                      Stage=pheno$`TCGA Subtype`))
head(imv.risk)


imv.risk1=imv.risk
imv.risk1=crbind2DataFrame(imv.risk1)
table(imv.risk1$Stage)
imv.risk1$Stage[imv.risk1$Stage=='I'|imv.risk1$Stage=='II']='I+II'
imv.risk1$Stage[imv.risk1$Stage=='III'|imv.risk1$Stage=='IV']='III+IV'

fig8b=ggsurvplot(fit=survfit( Surv(OS.time/12, OS) ~ Risktype,
                              data = imv.risk1[which(imv.risk1$Stage=='I+II'),]),
                 data=imv.risk1[which(imv.risk1$Stage=='I+II'),],
                 #surv.median.line = "hv",
                 conf.int = T,pval = T,fun = "pct",risk.table =T, palette = risktype.col,
                 title='IMvigor210 Stage I+II',
                 linetype = c("solid", "dashed","strata")[1],
                 #legend = c('top', 'bottom', 'left', 'right', 'none')[5],
                 legend = c(0.8,0.75), 
                 legend.title = "",
                 legend.labs = c("High","Low"))
fig8b
fig8b=mg_merge_plot(fig8b$plot,fig8b$table,nrow=2,heights = c(3,1),align = 'v')
fig8c=ggsurvplot(fit=survfit( Surv(OS.time/12, OS) ~ Risktype,
                              data = imv.risk1[which(imv.risk1$Stage=='III+IV'),]),
                 data=imv.risk1[which(imv.risk1$Stage=='III+IV'),],
                 conf.int =T,pval = T,fun = "pct",risk.table = T,palette = risktype.col,
                 title='IMvigor210 Stage III+IV',
                 linetype = c("solid", "dashed","strata")[1],
                 #legend = c('top', 'bottom', 'left', 'right', 'none')[5],
                 legend = c(0.8,0.75), 
                 legend.title = "",
                 legend.labs = c("High","Low"))
fig8c
fig8c=mg_merge_plot(fig8c$plot,fig8c$table,nrow=2,heights = c(3,1),align = 'v')
table(imv.risk$binaryResponse)
fig8d=imv.risk[which(imv.risk$binaryResponse!='NA'),] %>%
  ggplot(aes(x=binaryResponse, y=Riskscore,fill = binaryResponse)) +
  #geom_violin()+  
  scale_fill_manual(values = pal_nejm()(10)[3:4])+
  geom_boxplot()+
  theme_classic(base_size = 20)+
  ggpubr::stat_compare_means(aes(group=binaryResponse), label = "p.format", method = 'wilcox.test')+
  theme_classic()+
  theme(legend.position = 'none',axis.text = element_text(color = 'black'),
        axis.text.x = element_text(size = 12),axis.title.x = element_text(size = 15),
        axis.text.y = element_text(size = 12),axis.title.y = element_text(size = 15))
fig8d



fig8=mg_merge_plot(mg_merge_plot(fig8a,fig8b,ncol=2,labels = LETTERS[1:2]),
                   mg_merge_plot(fig8d,fig8c,ncol=2,labels = LETTERS[3:4]),
                   
                   nrow=2,ncol=1)
fig8
ggsave('08_risktype.treatment/fig8.pdf',height =10 ,width = 8)

