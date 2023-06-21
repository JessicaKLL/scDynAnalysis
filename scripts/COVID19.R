library(scDynAnalysis)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(pheatmap)
library(ggformula)

DATA<-readRDS("./PBMC_vaccine_CITE.rds")
load("./Splitted_covid19_xdonor.RData")

#Random Forest

RF_DA<-list()
for (i in 1:4) {
  RF_DA[[i]]<-RandomForest(as.factor(cell_type)~.,data=DA[[i]][,1:177],strata = "cell_type",ntree=200,split=T)
}
RF_DB<-list()
for (i in 1:4) {
  RF_DB[[i]]<-RandomForest(as.factor(cell_type)~.,data=DB[[i]][,1:177],strata = "cell_type",ntree=200,split=T)
}
RF_DE<-list()
for (i in 1:4) {
  RF_DE[[i]]<-RandomForest(as.factor(cell_type)~.,data=DE[[i]][,1:177],strata = "cell_type",ntree=200,split=T)
}
RF_DF<-list()
for (i in 1:4) {
  RF_DF[[i]]<-RandomForest(as.factor(cell_type)~.,data=DF[[i]][,1:177],strata = "cell_type",ntree=200,split=T)
}
RF_DJ<-list()
for (i in 1:4) {
  RF_DJ[[i]]<-RandomForest(as.factor(cell_type)~.,data=DJ[[i]][,1:177],strata = "cell_type",ntree=200,split=T)
}
RF_DM<-list()
for (i in 1:4) {
  RF_DM[[i]]<-RandomForest(as.factor(cell_type)~.,data=DM[[i]][,1:177],strata = "cell_type",ntree=200,split=T)
}

#Quality check

DA_quality<-list()
for (i in 1:length(RF_DA)) {
  DA_quality[[i]]<-model_quality(RF_DA[[i]]$RF_model,RF_DA[[i]]$Test_Set,RF_DA[[i]]$Test_Set$cell_type)
}
names(DA_quality)<-names(DA)
DB_quality<-list()
for (i in 1:length(RF_DB)) {
  DB_quality[[i]]<-model_quality(RF_DB[[i]]$RF_model,RF_DB[[i]]$Test_Set,RF_DB[[i]]$Test_Set$cell_type)
}
names(DB_quality)<-names(DB)
DE_quality<-list()
for (i in 1:length(RF_DE)) {
  DE_quality[[i]]<-model_quality(RF_DE[[i]]$RF_model,RF_DE[[i]]$Test_Set,RF_DE[[i]]$Test_Set$cell_type)
}
names(DE_quality)<-names(DE)
DF_quality<-list()
for (i in 1:length(RF_DF)) {
  DF_quality[[i]]<-model_quality(RF_DF[[i]]$RF_model,RF_DF[[i]]$Test_Set,RF_DF[[i]]$Test_Set$cell_type)
}
names(DF_quality)<-names(DF)
DJ_quality<-list()
for (i in 1:length(RF_DJ)) {
  DJ_quality[[i]]<-model_quality(RF_DJ[[i]]$RF_model,RF_DJ[[i]]$Test_Set,RF_DJ[[i]]$Test_Set$cell_type)
}
names(DJ_quality)<-names(DJ)
DM_quality<-list()
for (i in 1:length(RF_DM)) {
  DM_quality[[i]]<-model_quality(RF_DM[[i]]$RF_model,RF_DM[[i]]$Test_Set,RF_DM[[i]]$Test_Set$cell_type)
}
names(DM_quality)<-names(DM)

# Selection of relevant features

DA_impFeat<-list()
for (i in 1:length(DA)) {
  DA_impFeat[[i]]<-important_feat(RF_DA[[i]]$RF_model,data=DA[[i]],n=30, main=paste0(names(DA)[i],": Important Features"))
}
names(DA_impFeat)<-names(DA)
DA_selected_feat<-c()
for (x in 1:length(DA_impFeat)) {
  for (y in 1:length(DA_impFeat[[x]]$Selected)) {
    if (DA_impFeat[[x]]$Selected[y] %in% DA_selected_feat){
      next
    }
    else{
      DA_selected_feat<-append(DA_selected_feat,DA_impFeat[[x]]$Selected[y])
    }
  }
}
DB_impFeat<-list()
for (i in 1:length(DB)) {
  DB_impFeat[[i]]<-important_feat(RF_DB[[i]]$RF_model,data=DB[[i]],n=30, main=paste0(names(DB)[i],": Important Features"))
}
names(DB_impFeat)<-names(DB)
DB_selected_feat<-c()
for (x in 1:length(DB_impFeat)) {
  for (y in 1:length(DB_impFeat[[x]]$Selected)) {
    if (DB_impFeat[[x]]$Selected[y] %in% DB_selected_feat){
      next
    }
    else{
      DB_selected_feat<-append(DB_selected_feat,DB_impFeat[[x]]$Selected[y])
    }
  }
}
DE_impFeat<-list()
for (i in 1:length(DE)) {
  DE_impFeat[[i]]<-important_feat(RF_DE[[i]]$RF_model,data=DE[[i]],n=30, main=paste0(names(DE)[i],": Important Features"))
}
names(DE_impFeat)<-names(DE)
DE_selected_feat<-c()
for (x in 1:length(DE_impFeat)) {
  for (y in 1:length(DE_impFeat[[x]]$Selected)) {
    if (DE_impFeat[[x]]$Selected[y] %in% DE_selected_feat){
      next
    }
    else{
      DE_selected_feat<-append(DE_selected_feat,DE_impFeat[[x]]$Selected[y])
    }
  }
}
DF_impFeat<-list()
for (i in 1:length(DF)) {
  DF_impFeat[[i]]<-important_feat(RF_DF[[i]]$RF_model,data=DF[[i]],n=30, main=paste0(names(DF)[i],": Important Features"))
}
names(DF_impFeat)<-names(DF)
DF_selected_feat<-c()
for (x in 1:length(DF_impFeat)) {
  for (y in 1:length(DF_impFeat[[x]]$Selected)) {
    if (DF_impFeat[[x]]$Selected[y] %in% DF_selected_feat){
      next
    }
    else{
      DF_selected_feat<-append(DF_selected_feat,DF_impFeat[[x]]$Selected[y])
    }
  }
}
DJ_impFeat<-list()
for (i in 1:length(DJ)) {
  DJ_impFeat[[i]]<-important_feat(RF_DJ[[i]]$RF_model,data=DJ[[i]],n=30, main=paste0(names(DJ)[i],": Important Features"))
}
names(DJ_impFeat)<-names(DJ)
DJ_selected_feat<-c()
for (x in 1:length(DJ_impFeat)) {
  for (y in 1:length(DJ_impFeat[[x]]$Selected)) {
    if (DJ_impFeat[[x]]$Selected[y] %in% DJ_selected_feat){
      next
    }
    else{
      DJ_selected_feat<-append(DJ_selected_feat,DJ_impFeat[[x]]$Selected[y])
    }
  }
}

DM_impFeat<-list()
for (i in 1:length(DM)) {
  DM_impFeat[[i]]<-important_feat(RF_DM[[i]]$RF_model,data=DM[[i]],n=30, main=paste0(names(DM)[i],": Important Features"))
}
names(DM_impFeat)<-names(DM)
DM_selected_feat<-c()
for (x in 1:length(DM_impFeat)) {
  for (y in 1:length(DM_impFeat[[x]]$Selected)) {
    if (DM_impFeat[[x]]$Selected[y] %in% DM_selected_feat){
      next
    }
    else{
      DM_selected_feat<-append(DM_selected_feat,DM_impFeat[[x]]$Selected[y])
    }
  }
}

# First filtering

for (i in 1:length(DA)) {
  cell_type<-DA[[i]]$cell_type
  DA[[i]]<-DA[[i]][,which(colnames(DA[[i]]) %in% Selected_feat)]
  DA[[i]]$cell_type<-cell_type
}
for (i in 1:length(DB)) {
  cell_type<-DB[[i]]$cell_type
  DB[[i]]<-DB[[i]][,which(colnames(DB[[i]]) %in% Selected_feat)]
  DB[[i]]$cell_type<-cell_type
}
for (i in 1:length(DE)) {
  cell_type<-DE[[i]]$cell_type
  DE[[i]]<-DE[[i]][,which(colnames(DE[[i]]) %in% Selected_feat)]
  DE[[i]]$cell_type<-cell_type
}
for (i in 1:length(DF)) {
  cell_type<-DF[[i]]$cell_type
  DF[[i]]<-DF[[i]][,which(colnames(DF[[i]]) %in% Selected_feat)]
  DF[[i]]$cell_type<-cell_type
}
for (i in 1:length(DJ)) {
  cell_type<-DJ[[i]]$cell_type
  DJ[[i]]<-DJ[[i]][,which(colnames(DJ[[i]]) %in% Selected_feat)]
  DJ[[i]]$cell_type<-cell_type
}
for (i in 1:length(DM)) {
  cell_type<-DM[[i]]$cell_type
  DM[[i]]<-DM[[i]][,which(colnames(DM[[i]]) %in% Selected_feat)]
  DM[[i]]$cell_type<-cell_type
}

# Technical variance signicance

suppressWarnings(DA_Within_var_0_2<-Within_variance(DA$Day0,DA$Day2,n=100,perc=0.3))
suppressWarnings(DA_Within_var_2_11<-Within_variance(DA$Day2,DA$Day11,n=100,perc=0.3))
suppressWarnings(DA_Within_var_11_28<-Within_variance(DA$Day11,DA$Day28,n=100,perc=0.3))
DA_p_values<-data.frame(DA_0_2=mean(DA_Within_var_0_2$p.value),
                        DA_2_11=mean(DA_Within_var_2_11$p.value),
                        DA_11_28=mean(DA_Within_var_11_28$p.value))

suppressWarnings(DB_Within_var_0_2<-Within_variance(DB$Day0,DB$Day2,n=100,perc=0.3))
suppressWarnings(DB_Within_var_2_11<-Within_variance(DB$Day2,DB$Day11,n=100,perc=0.3))
suppressWarnings(DB_Within_var_11_28<-Within_variance(DB$Day11,DB$Day28,n=100,perc=0.3))
DB_p_values<-data.frame(DB_0_2=mean(DB_Within_var_0_2$p.value),
                        DB_2_11=mean(DB_Within_var_2_11$p.value),
                        DB_11_28=mean(DB_Within_var_11_28$p.value))

suppressWarnings(DE_Within_var_0_2<-Within_variance(DE$Day0,DE$Day2,n=100,perc=0.3))
suppressWarnings(DE_Within_var_2_11<-Within_variance(DE$Day2,DE$Day11,n=100,perc=0.3))
suppressWarnings(DE_Within_var_11_28<-Within_variance(DE$Day11,DE$Day28,n=100,perc=0.3))
DE_p_values<-data.frame(DE_0_2=mean(DE_Within_var_0_2$p.value),
                        DE_2_11=mean(DE_Within_var_2_11$p.value),
                        DE_11_28=mean(DE_Within_var_11_28$p.value))

suppressWarnings(DF_Within_var_0_2<-Within_variance(DF$Day0,DF$Day2,n=100,perc=0.3))
suppressWarnings(DF_Within_var_2_11<-Within_variance(DF$Day2,DF$Day11,n=100,perc=0.3))
suppressWarnings(DF_Within_var_11_28<-Within_variance(DF$Day11,DF$Day28,n=100,perc=0.3))
DF_p_values<-data.frame(DF_0_2=mean(DF_Within_var_0_2$p.value),
                        DF_2_11=mean(DF_Within_var_2_11$p.value),
                        DF_11_28=mean(DF_Within_var_11_28$p.value))

suppressWarnings(DJ_Within_var_0_2<-Within_variance(DJ$Day0,DJ$Day2,n=100,perc=0.3))
suppressWarnings(DJ_Within_var_2_11<-Within_variance(DJ$Day2,DJ$Day11,n=100,perc=0.3))
suppressWarnings(DJ_Within_var_11_28<-Within_variance(DJ$Day11,DJ$Day28,n=100,perc=0.3))
DJ_p_values<-data.frame(DJ_0_2=mean(DJ_Within_var_0_2$p.value),
                        DJ_2_11=mean(DJ_Within_var_2_11$p.value),
                        DJ_11_28=mean(DJ_Within_var_11_28$p.value))

suppressWarnings(DM_Within_var_0_2<-Within_variance(DM$Day0,DM$Day2,n=100,perc=0.3))
suppressWarnings(DM_Within_var_2_11<-Within_variance(DM$Day2,DM$Day11,n=100,perc=0.3))
suppressWarnings(DM_Within_var_11_28<-Within_variance(DM$Day11,DM$Day28,n=100,perc=0.3))
DM_p_values<-data.frame(DM_0_2=mean(DM_Within_var_0_2$p.value),
                        DM_2_11=mean(DM_Within_var_2_11$p.value),
                        DM_11_28=mean(DM_Within_var_11_28$p.value))

# Altered features

DA_0<-DA$Day0
DA_2<-DA$Day2
DA_11<-DA$Day11
DA_28<-DA$Day28
DA_0$time_point<-"Day0"
DA_2$time_point<-"Day2"
DA_11$time_point<-"Day11"
DA_28$time_point<-"Day28"
DA_complete<-rbind(DA_0,DA_2)
DA_complete<-rbind(DA_complete,DA_11)
DA_complete<-rbind(DA_complete,DA_28)
DA_Between_var<-Between_variance(DA_complete,Selected_feat,DA_complete$time_point)
DA_pVal_alter<-select_alter(DA_Between_var,Selected_feat,explicit = T)

DB_0<-DB$Day0
DB_2<-DB$Day2
DB_11<-DB$Day11
DB_28<-DB$Day28
DB_0$time_point<-"Day0"
DB_2$time_point<-"Day2"
DB_11$time_point<-"Day11"
DB_28$time_point<-"Day28"
DB_complete<-rbind(DB_0,DB_2)
DB_complete<-rbind(DB_complete,DB_11)
DB_complete<-rbind(DB_complete,DB_28)
DB_Between_var<-Between_variance(DB_complete,Selected_feat,DB_complete$time_point)
DB_pVal_alter<-select_alter(DB_Between_var,Selected_feat,explicit = T)

DE_0<-DE$Day0
DE_2<-DE$Day2
DE_11<-DE$Day11
DE_28<-DE$Day28
DE_0$time_point<-"Day0"
DE_2$time_point<-"Day2"
DE_11$time_point<-"Day11"
DE_28$time_point<-"Day28"
DE_complete<-rbind(DE_0,DE_2)
DE_complete<-rbind(DE_complete,DE_11)
DE_complete<-rbind(DE_complete,DE_28)
DE_Between_var<-Between_variance(DE_complete,Selected_feat,DE_complete$time_point)
DE_pVal_alter<-select_alter(DE_Between_var,Selected_feat,explicit = T)

DF_0<-DF$Day0
DF_2<-DF$Day2
DF_11<-DF$Day11
DF_28<-DF$Day28
DF_0$time_point<-"Day0"
DF_2$time_point<-"Day2"
DF_11$time_point<-"Day11"
DF_28$time_point<-"Day28"
DF_complete<-rbind(DF_0,DF_2)
DF_complete<-rbind(DF_complete,DF_11)
DF_complete<-rbind(DF_complete,DF_28)
DF_Between_var<-Between_variance(DF_complete,Selected_feat,DF_complete$time_point)
DF_pVal_alter<-select_alter(DF_Between_var,Selected_feat,explicit = T)

DJ_0<-DJ$Day0
DJ_2<-DJ$Day2
DJ_11<-DJ$Day11
DJ_28<-DJ$Day28
DJ_0$time_point<-"Day0"
DJ_2$time_point<-"Day2"
DJ_11$time_point<-"Day11"
DJ_28$time_point<-"Day28"
DJ_complete<-rbind(DJ_0,DJ_2)
DJ_complete<-rbind(DJ_complete,DJ_11)
DJ_complete<-rbind(DJ_complete,DJ_28)
DJ_Between_var<-Between_variance(DJ_complete,Selected_feat,DJ_complete$time_point)
DJ_pVal_alter<-select_alter(DJ_Between_var,Selected_feat,explicit = T)

DM_0<-DM$Day0
DM_2<-DM$Day2
DM_11<-DM$Day11
DM_28<-DM$Day28
DM_0$time_point<-"Day0"
DM_2$time_point<-"Day2"
DM_11$time_point<-"Day11"
DM_28$time_point<-"Day28"
DM_complete<-rbind(DM_0,DM_2)
DM_complete<-rbind(DM_complete,DM_11)
DM_complete<-rbind(DM_complete,DM_28)
DM_Between_var<-Between_variance(DM_complete,Selected_feat,DM_complete$time_point)
DM_pVal_alter<-select_alter(DM_Between_var,Selected_feat,explicit = T)

pVal_alter<-c(DA_pVal_alter,DB_pVal_alter,DE_pVal_alter,DF_pVal_alter,DJ_pVal_alter,DM_pVal_alter,"antigen")
pVal_alter<-unique(pVal_alter)

# Second filtering

for (i in 1:length(DA)) {
  cell_type<-DA[[i]]$cell_type
  DA[[i]]<-DA[[i]][,which(colnames(DA[[i]]) %in% pVal_alter)]
  DA[[i]]$cell_type<-cell_type
}
for (i in 1:length(DB)) {
  cell_type<-DB[[i]]$cell_type
  DB[[i]]<-DB[[i]][,which(colnames(DB[[i]]) %in% pVal_alter)]
  DB[[i]]$cell_type<-cell_type
}
for (i in 1:length(DE)) {
  cell_type<-DE[[i]]$cell_type
  DE[[i]]<-DE[[i]][,which(colnames(DE[[i]]) %in% pVal_alter)]
  DE[[i]]$cell_type<-cell_type
}
for (i in 1:length(DF)) {
  cell_type<-DF[[i]]$cell_type
  DF[[i]]<-DF[[i]][,which(colnames(DF[[i]]) %in% pVal_alter)]
  DF[[i]]$cell_type<-cell_type
}
for (i in 1:length(DJ)) {
  cell_type<-DJ[[i]]$cell_type
  DJ[[i]]<-DJ[[i]][,which(colnames(DJ[[i]]) %in% pVal_alter)]
  DJ[[i]]$cell_type<-cell_type
}
for (i in 1:length(DM)) {
  cell_type<-DM[[i]]$cell_type
  DM[[i]]<-DM[[i]][,which(colnames(DM[[i]]) %in% pVal_alter)]
  DM[[i]]$cell_type<-cell_type
}

## feature plots
Features<-pVal_alter[-which(pVal_alter=="antigen")]
P<-FeaturePlot(DATA, features = Features,cols = c("lightgrey","red2"),split.by = "timepoint")

## boxplots

tp<-list()
tp[[1]]<-compare_tp(DA_complete,Cell_Type = DA_complete$cell_type,Time_Step = DA_complete$time_point,Feature="CD38-1",main="")
tp[[2]]<-compare_tp(DB_complete,Cell_Type = DB_complete$cell_type,Time_Step = DB_complete$time_point,Feature="CD38-1",main="")
tp[[3]]<-compare_tp(DE_complete,Cell_Type = DE_complete$cell_type,Time_Step = DE_complete$time_point,Feature="CD38-1",main="")
tp[[4]]<-compare_tp(DF_complete,Cell_Type = DF_complete$cell_type,Time_Step = DF_complete$time_point,Feature="CD38-1",main="")
tp[[5]]<-compare_tp(DJ_complete,Cell_Type = DJ_complete$cell_type,Time_Step = DJ_complete$time_point,Feature="CD38-1",main="")
tp[[6]]<-compare_tp(DM_complete,Cell_Type = DM_complete$cell_type,Time_Step = DM_complete$time_point,Feature="CD38-1",main="")
CD38<-ggarrange(tp[[1]],tp[[2]],tp[[3]],tp[[4]],tp[[5]],tp[[6]],
                nrow = 1,ncol = 6,common.legend = T,legend = "right",labels = c("DA","DB","DE","DF","DJ","DM"))
CD38<-annotate_figure(CD38, top = text_grob("CD38-1",color = "black", face = "bold", size = 14))

tp<-list()
tp[[1]]<-compare_tp(P1_complete,Cell_Type = P1_complete$cell_type,Time_Step = P1_complete$time_point,Feature="GZMB",main="")
tp[[2]]<-compare_tp(P2_complete,Cell_Type = P2_complete$cell_type,Time_Step = P2_complete$time_point,Feature="GZMB",main="")
tp[[3]]<-compare_tp(P3_complete,Cell_Type = P3_complete$cell_type,Time_Step = P3_complete$time_point,Feature="GZMB",main="")
tp[[4]]<-compare_tp(P4_complete,Cell_Type = P4_complete$cell_type,Time_Step = P4_complete$time_point,Feature="GZMB",main="")
tp[[5]]<-compare_tp(P5_complete,Cell_Type = P5_complete$cell_type,Time_Step = P5_complete$time_point,Feature="GZMB",main="")
GZMB<-ggarrange(plotlist = tp,nrow = 1,ncol = 5,common.legend = T,legend = "right",labels = c("P1","P2","P3","P4","P5"))
GZMB<-annotate_figure(GZMB, top = text_grob("GZMB",color = "black", face = "bold", size = 14))

tp<-list()
tp[[1]]<-compare_tp(P1_complete,Cell_Type = P1_complete$cell_type,Time_Step = P1_complete$time_point,Feature="TRDV1",main="")
tp[[2]]<-compare_tp(P2_complete,Cell_Type = P2_complete$cell_type,Time_Step = P2_complete$time_point,Feature="TRDV1",main="")
tp[[3]]<-compare_tp(P3_complete,Cell_Type = P3_complete$cell_type,Time_Step = P3_complete$time_point,Feature="TRDV1",main="")
tp[[4]]<-compare_tp(P4_complete,Cell_Type = P4_complete$cell_type,Time_Step = P4_complete$time_point,Feature="TRDV1",main="")
tp[[5]]<-compare_tp(P5_complete,Cell_Type = P5_complete$cell_type,Time_Step = P5_complete$time_point,Feature="TRDV1",main="")
TRDV1<-ggarrange(plotlist = tp,nrow = 1,ncol = 5,common.legend = T,legend = "right",labels = c("P1","P2","P3","P4","P5"))
TRDV1<-annotate_figure(TRDV1, top = text_grob("TRDV1",color = "black", face = "bold", size = 14))

tp<-list()
tp[[1]]<-compare_tp(P1_complete,Cell_Type = P1_complete$cell_type,Time_Step = P1_complete$time_point,Feature="EOMES",main="")
tp[[2]]<-compare_tp(P2_complete,Cell_Type = P2_complete$cell_type,Time_Step = P2_complete$time_point,Feature="EOMES",main="")
tp[[3]]<-compare_tp(P3_complete,Cell_Type = P3_complete$cell_type,Time_Step = P3_complete$time_point,Feature="EOMES",main="")
tp[[4]]<-compare_tp(P4_complete,Cell_Type = P4_complete$cell_type,Time_Step = P4_complete$time_point,Feature="EOMES",main="")
tp[[5]]<-compare_tp(P5_complete,Cell_Type = P5_complete$cell_type,Time_Step = P5_complete$time_point,Feature="EOMES",main="")
EOMES<-ggarrange(plotlist = tp,nrow = 1,ncol = 5,common.legend = T,legend = "right",labels = c("P1","P2","P3","P4","P5"))
EOMES<-annotate_figure(EOMES, top = text_grob("EOMES",color = "black", face = "bold", size = 14))

tp<-list()
tp[[1]]<-compare_tp(P1_complete,Cell_Type = P1_complete$cell_type,Time_Step = P1_complete$time_point,Feature="BATF",main="")
tp[[2]]<-compare_tp(P2_complete,Cell_Type = P2_complete$cell_type,Time_Step = P2_complete$time_point,Feature="BATF",main="")
tp[[3]]<-compare_tp(P3_complete,Cell_Type = P3_complete$cell_type,Time_Step = P3_complete$time_point,Feature="BATF",main="")
tp[[4]]<-compare_tp(P4_complete,Cell_Type = P4_complete$cell_type,Time_Step = P4_complete$time_point,Feature="BATF",main="")
tp[[5]]<-compare_tp(P5_complete,Cell_Type = P5_complete$cell_type,Time_Step = P5_complete$time_point,Feature="BATF",main="")
BATF<-ggarrange(plotlist = tp,nrow = 1,ncol = 5,common.legend = T,legend = "right",labels = c("P1","P2","P3","P4","P5"))
BATF<-annotate_figure(BATF, top = text_grob("BATF",color = "black", face = "bold", size = 14))

# cell state ordering

DA_0_cell_type<-DA$Day0$cell_type
DA_2_cell_type<-DA$Day2$cell_type
DA_11_cell_type<-DA$Day11$cell_type
DA_28_cell_type<-DA$Day28$cell_type

DB_0_cell_type<-DB$Day0$cell_type
DB_2_cell_type<-DB$Day2$cell_type
DB_11_cell_type<-DB$Day11$cell_type
DB_28_cell_type<-DB$Day28$cell_type

DE_0_cell_type<-DE$Day0$cell_type
DE_2_cell_type<-DE$Day2$cell_type
DE_11_cell_type<-DE$Day11$cell_type
DE_28_cell_type<-DE$Day28$cell_type

DF_0_cell_type<-DF$Day0$cell_type
DF_2_cell_type<-DF$Day2$cell_type
DF_11_cell_type<-DF$Day11$cell_type
DF_28_cell_type<-DF$Day28$cell_type

DJ_0_cell_type<-DJ$Day0$cell_type
DJ_2_cell_type<-DJ$Day2$cell_type
DJ_11_cell_type<-DJ$Day11$cell_type
DJ_28_cell_type<-DJ$Day28$cell_type

DM_0_cell_type<-DM$Day0$cell_type
DM_2_cell_type<-DM$Day2$cell_type
DM_11_cell_type<-DM$Day11$cell_type
DM_28_cell_type<-DM$Day28$cell_type

DA$Day0$cell_type<-NULL
DA$Day2$cell_type<-NULL
DA$Day11$cell_type<-NULL
DA$Day28$cell_type<-NULL

DB$Day0$cell_type<-NULL
DB$Day2$cell_type<-NULL
DB$Day11$cell_type<-NULL
DB$Day28$cell_type<-NULL

DE$Day0$cell_type<-NULL
DE$Day2$cell_type<-NULL
DE$Day11$cell_type<-NULL
DE$Day28$cell_type<-NULL

DF$Day0$cell_type<-NULL
DF$Day2$cell_type<-NULL
DF$Day11$cell_type<-NULL
DF$Day28$cell_type<-NULL

DJ$Day0$cell_type<-NULL
DJ$Day2$cell_type<-NULL
DJ$Day11$cell_type<-NULL
DJ$Day28$cell_type<-NULL

DM$Day0$cell_type<-NULL
DM$Day2$cell_type<-NULL
DM$Day11$cell_type<-NULL
DM$Day28$cell_type<-NULL

DA$Day0$time_point<-NULL
DA$Day2$time_point<-NULL
DA$Day11$time_point<-NULL
DA$Day28$time_point<-NULL

DB$Day0$time_point<-NULL
DB$Day2$time_point<-NULL
DB$Day11$time_point<-NULL
DB$Day28$time_point<-NULL

DE$Day0$time_point<-NULL
DE$Day2$time_point<-NULL
DE$Day11$time_point<-NULL
DE$Day28$time_point<-NULL

DF$Day0$time_point<-NULL
DF$Day2$time_point<-NULL
DF$Day11$time_point<-NULL
DF$Day28$time_point<-NULL

DJ$Day0$time_point<-NULL
DJ$Day2$time_point<-NULL
DJ$Day11$time_point<-NULL
DJ$Day28$time_point<-NULL

DM$Day0$time_point<-NULL
DM$Day2$time_point<-NULL
DM$Day11$time_point<-NULL
DM$Day28$time_point<-NULL

DA_hclust<-divHier(DA,k=100)
DB_hclust<-divHier(DB,k=100)
DE_hclust<-divHier(DE,k=100)
DF_hclust<-divHier(DF,k=100)
DJ_hclust<-divHier(DJ,k=100)
DM_hclust<-divHier(DM,k=100)

DA_new_df<-gen_new_data(DA_hclust,method = "sum")
DB_new_df<-gen_new_data(DB_hclust,method = "sum")
DE_new_df<-gen_new_data(DE_hclust,method = "sum")
DF_new_df<-gen_new_data(DF_hclust,method = "sum")
DJ_new_df<-gen_new_data(DJ_hclust,method = "sum")
DM_new_df<-gen_new_data(DM_hclust,method="sum")

feat_num<-length(pVal_alter)

DA<-split(DA_new_df,DA_new_df$time_point)
DA0_2<-rbind(DA$Day0,DA$Day2)
DA2_11<-rbind(DA$Day2,DA$Day11)
DA11_28<-rbind(DA$Day11,DA$Day28)

DA_corr0_2<-find_transition(DA0_2,cluster = "clusters",feat_num = feat_num)
DA_corr2_11<-find_transition(DA2_11,cluster = "clusters",feat_num = feat_num)
DA_corr11_28<-find_transition(DA11_28,cluster = "clusters",feat_num = feat_num)

DB<-split(DB_new_df,DB_new_df$time_point)
DB0_2<-rbind(DB$Day0,DB$Day2)
DB2_11<-rbind(DB$Day2,DB$Day11)
DB11_28<-rbind(DB$Day11,DB$Day28)

DB_corr0_2<-find_transition(DB0_2,cluster = "clusters",feat_num = feat_num)
DB_corr2_11<-find_transition(DB2_11,cluster = "clusters",feat_num = feat_num)
DB_corr11_28<-find_transition(DB11_28,cluster = "clusters",feat_num = feat_num)

DE<-split(DE_new_df,DE_new_df$time_point)
DE0_2<-rbind(DE$Day0,DE$Day2)
DE2_11<-rbind(DE$Day2,DE$Day11)
DE11_28<-rbind(DE$Day11,DE$Day28)

DE_corr0_2<-find_transition(DE0_2,cluster = "clusters",feat_num = feat_num)
DE_corr2_11<-find_transition(DE2_11,cluster = "clusters",feat_num = feat_num)
DE_corr11_28<-find_transition(DE11_28,cluster = "clusters",feat_num = feat_num)

DF<-split(DF_new_df,DF_new_df$time_point)
DF0_2<-rbind(DF$Day0,DF$Day2)
DF2_11<-rbind(DF$Day2,DF$Day11)
DF11_28<-rbind(DF$Day11,DF$Day28)

DF_corr0_2<-find_transition(DF0_2,cluster = "clusters",feat_num = feat_num)
DF_corr2_11<-find_transition(DF2_11,cluster = "clusters",feat_num = feat_num)
DF_corr11_28<-find_transition(DF11_28,cluster = "clusters",feat_num = feat_num)

DJ<-split(DJ_new_df,DJ_new_df$time_point)
DJ0_2<-rbind(DJ$Day0,DJ$Day2)
DJ2_11<-rbind(DJ$Day2,DJ$Day11)
DJ11_28<-rbind(DJ$Day11,DJ$Day28)

DJ_corr0_2<-find_transition(DJ0_2,cluster = "clusters",feat_num = feat_num)
DJ_corr2_11<-find_transition(DJ2_11,cluster = "clusters",feat_num = feat_num)
DJ_corr11_28<-find_transition(DJ11_28,cluster = "clusters",feat_num = feat_num)

DM<-split(DM_new_df,DM_new_df$time_point)
DM0_2<-rbind(DM$Day0,DM$Day2)
DM2_11<-rbind(DM$Day2,DM$Day11)
DM11_28<-rbind(DM$Day11,DM$Day28)

DM_corr0_2<-find_transition(DM0_2,cluster = "clusters",feat_num = feat_num)
DM_corr2_11<-find_transition(DM2_11,cluster = "clusters",feat_num = feat_num)
DM_corr11_28<-find_transition(DM11_28,cluster = "clusters",feat_num = feat_num)

DA_0_order<-cell_state_order(DA_new_df[which(DA_new_df$time_point=="Day0"),],root = DA_corr0_2$Identical$x,feat_num = feat_num,cluster = "clusters")
DA_2_order<-cell_state_order(DA_new_df[which(DA_new_df$time_point=="Day2"),],root = DA_corr0_2$Identical$y,feat_num = feat_num,cluster = "clusters")
DA_11_order<-cell_state_order(DA_new_df[which(DA_new_df$time_point=="Day11"),],root = DA_corr11_28$Identical$x,feat_num = feat_num,cluster = "clusters")
DA_28_order<-cell_state_order(DA_new_df[which(DA_new_df$time_point=="Day28"),],root = DA_corr11_28$Identical$y,feat_num = feat_num,cluster = "clusters")
DA_order<-c(rev(DA_0_order),DA_2_order,rev(DA_11_order),DA_28_order)

DB_0_order<-cell_state_order(DB_new_df[which(DB_new_df$time_point=="Day0"),],root = DB_corr0_2$Identical$x,feat_num = feat_num,cluster = "clusters")
DB_2_order<-cell_state_order(DB_new_df[which(DB_new_df$time_point=="Day2"),],root = DB_corr0_2$Identical$y,feat_num = feat_num,cluster = "clusters")
DB_11_order<-cell_state_order(DB_new_df[which(DB_new_df$time_point=="Day11"),],root = DB_corr11_28$Identical$x,feat_num = feat_num,cluster = "clusters")
DB_28_order<-cell_state_order(DB_new_df[which(DB_new_df$time_point=="Day28"),],root = DB_corr11_28$Identical$y,feat_num = feat_num,cluster = "clusters")
DB_order<-c(rev(DB_0_order),DB_2_order,rev(DB_11_order),DB_28_order)

DE_0_order<-cell_state_order(DE_new_df[which(DE_new_df$time_point=="Day0"),],root = DE_corr0_2$Identical$x,feat_num = feat_num,cluster = "clusters")
DE_2_order<-cell_state_order(DE_new_df[which(DE_new_df$time_point=="Day2"),],root = DE_corr0_2$Identical$y,feat_num = feat_num,cluster = "clusters")
DE_11_order<-cell_state_order(DE_new_df[which(DE_new_df$time_point=="Day11"),],root = DE_corr11_28$Identical$x,feat_num = feat_num,cluster = "clusters")
DE_28_order<-cell_state_order(DE_new_df[which(DE_new_df$time_point=="Day28"),],root = DE_corr11_28$Identical$y,feat_num = feat_num,cluster = "clusters")
DE_order<-c(rev(DE_0_order),DE_2_order,rev(DE_11_order),DE_28_order)

DF_0_order<-cell_state_order(DF_new_df[which(DF_new_df$time_point=="Day0"),],root = DF_corr0_2$Identical$x,feat_num = feat_num,cluster = "clusters")
DF_2_order<-cell_state_order(DF_new_df[which(DF_new_df$time_point=="Day2"),],root = DF_corr0_2$Identical$y,feat_num = feat_num,cluster = "clusters")
DF_11_order<-cell_state_order(DF_new_df[which(DF_new_df$time_point=="Day11"),],root = DF_corr11_28$Identical$x,feat_num = feat_num,cluster = "clusters")
DF_28_order<-cell_state_order(DF_new_df[which(DF_new_df$time_point=="Day28"),],root = DF_corr11_28$Identical$y,feat_num = feat_num,cluster = "clusters")
DF_order<-c(rev(DF_0_order),DF_2_order,rev(DF_11_order),DF_28_order)

DJ_0_order<-cell_state_order(DJ_new_df[which(DJ_new_df$time_point=="Day0"),],root = DJ_corr0_2$Identical$x,feat_num = feat_num,cluster = "clusters")
DJ_2_order<-cell_state_order(DJ_new_df[which(DJ_new_df$time_point=="Day2"),],root = DJ_corr0_2$Identical$y,feat_num = feat_num,cluster = "clusters")
DJ_11_order<-cell_state_order(DJ_new_df[which(DJ_new_df$time_point=="Day11"),],root = DJ_corr11_28$Identical$x,feat_num = feat_num,cluster = "clusters")
DJ_28_order<-cell_state_order(DJ_new_df[which(DJ_new_df$time_point=="Day28"),],root = DJ_corr11_28$Identical$y,feat_num = feat_num,cluster = "clusters")
DJ_order<-c(rev(DJ_0_order),DJ_2_order,rev(DJ_11_order),DJ_28_order)

DM_0_order<-cell_state_order(DM_new_df[which(DM_new_df$time_point=="Day0"),],root = DM_corr0_2$Identical$x,feat_num = feat_num,cluster = "clusters")
DM_2_order<-cell_state_order(DM_new_df[which(DM_new_df$time_point=="Day2"),],root = DM_corr0_2$Identical$y,feat_num = feat_num,cluster = "clusters")
DM_11_order<-cell_state_order(DM_new_df[which(DM_new_df$time_point=="Day11"),],root = DM_corr11_28$Identical$x,feat_num = feat_num,cluster = "clusters")
DM_28_order<-cell_state_order(DM_new_df[which(DM_new_df$time_point=="Day28"),],root = DM_corr11_28$Identical$y,feat_num = feat_num,cluster = "clusters")
DM_order<-c(rev(DM_0_order),DM_2_order,rev(DM_11_order),DM_28_order)

DA_ordered_new_df<-DA_new_df[match(DA_order,DA_new_df$clusters),]
DB_ordered_new_df<-DB_new_df[match(DB_order,DB_new_df$clusters),]
DE_ordered_new_df<-DE_new_df[match(DE_order,DE_new_df$clusters),]
DF_ordered_new_df<-DF_new_df[match(DF_order,DF_new_df$clusters),]
DJ_ordered_new_df<-DJ_new_df[match(DJ_order,DJ_new_df$clusters),]
DM_ordered_new_df<-DM_new_df[match(DM_order,DM_new_df$clusters),]

# Dynamics

## spaghetti
ordered_donors_new_df<-list(DonorA=DA_ordered_new_df,DonorB=DB_ordered_new_df,DonorE=DE_ordered_new_df,
                            DonorF=DF_ordered_new_df,DonorJ=DJ_ordered_new_df,DonorM=DM_ordered_new_df)
ot_Gene_dyn_Xpatient<-list()
for (i in 1:length(Features)) {
  p1<-as.data.frame(ordered_donors_new_df[[1]][,Features[i]])
  p1$order<-ordered_donors_new_df[[1]]$clusters
  p1$donor<-"DonorA"
  names(p1)<-c("Gene","Meta_cells","Individual")
  p2<-as.data.frame(ordered_donors_new_df[[2]][,Features[i]])
  p2$order<-ordered_donors_new_df[[2]]$clusters
  p2$donor<-"DonorB"
  names(p2)<-c("Gene","Meta_cells","Individual")
  p3<-as.data.frame(ordered_donors_new_df[[3]][,Features[i]])
  p3$order<-ordered_donors_new_df[[3]]$clusters
  p3$donor<-"DonorE"
  names(p3)<-c("Gene","Meta_cells","Individual")
  p4<-as.data.frame(ordered_donors_new_df[[4]][,Features[i]])
  p4$order<-ordered_donors_new_df[[4]]$clusters
  p4$donor<-"DonorF"
  names(p4)<-c("Gene","Meta_cells","Individual")
  p5<-as.data.frame(ordered_donors_new_df[[5]][,Features[i]])
  p5$order<-ordered_donors_new_df[[5]]$clusters
  p5$donor<-"DonorJ"
  names(p5)<-c("Gene","Meta_cells","Individual")
  p6<-as.data.frame(ordered_donors_new_df[[6]][,Features[i]])
  p6$order<-ordered_donors_new_df[[6]]$clusters
  p6$donor<-"DonorM"
  names(p6)<-c("Gene","Meta_cells","Individual")
  df<-rbind(p1,p2)
  df<-rbind(df,p3)
  df<-rbind(df,p4)
  df<-rbind(df,p5)
  df<-rbind(df,p6)
  df$Meta_cells <- factor(df$Meta_cells, levels=unique(df$Meta_cells))
  df$Patient<-as.factor(df$donor)
  p<-plot_dyn(df)+
      theme_classic()+theme(axis.text.x = element_blank())+ylab(paste0(Features[i]))+
      labs(title = paste0(Features[i]))+
    geom_vline(xintercept = c(100), linetype=2, size = 0.3)+
    geom_vline(xintercept = c(200), linetype=2, size = 0.3)+
    geom_vline(xintercept = c(300), linetype=2, size = 0.3)+
    geom_text(aes(x=50, label="\nDay0", y=10), colour="black") +
    geom_text(aes(x=150, label="\nDay2", y=10), colour="black")+
    geom_text(aes(x=250, label="\nDay11", y=10), colour="black")+
    geom_text(aes(x=350, label="\nDay28", y=10), colour="black")+
    scale_color_manual(values = c(DonorA="red",DonorB="purple",DonorE="orange",
                                  DonorF="blue2",DonorJ="yellowgreen",DonorM="turquoise"))+
    theme(plot.title = element_text(hjust = 0.5))+ylab("")+xlab("")
  ot_Gene_dyn_Xpatient[[i]]<-p
}
names(ot_Gene_dyn_Xpatient)<-Features

## heatmaps

heatmaps<-list()
heatmaps[[1]]<-Expr_Heatmap(DA_ordered_new_df,Features,main = "DA")
heatmaps[[2]]<-Expr_Heatmap(DB_ordered_new_df,Features,main = "DB")
heatmaps[[3]]<-Expr_Heatmap(DE_ordered_new_df,Features,main = "DE")
heatmaps[[4]]<-Expr_Heatmap(DF_ordered_new_df,Features,main = "DF")
heatmaps[[5]]<-Expr_Heatmap(DJ_ordered_new_df,Features,main = "DJ")
heatmaps[[6]]<-Expr_Heatmap(DM_ordered_new_df,Features,main = "DM")

# Forecasting

Features<-pVal_alter[-which(pVal_alter=="antigen")]

pred_genes<-list()
for (i in 1:length(Features)) {
  output_1<-Cross_Validation_lstm(DA_ordered_new_df,Feature = Features[i],epochs = 1000) %>%
    plot_prediction(alpha = 0.65) +
    theme(legend.position = "bottom")
  output_2<-Cross_Validation_lstm(DB_ordered_new_df,Feature = Features[i],epochs = 1000) %>%
    plot_prediction(alpha = 0.65) +
    theme(legend.position = "bottom")
  output_3<-Cross_Validation_lstm(DE_ordered_new_df,Feature = Features[i],epochs = 1000) %>%
    plot_prediction(alpha = 0.65) +
    theme(legend.position = "bottom")
  output_4<-Cross_Validation_lstm(DF_ordered_new_df,Feature = Features[i],epochs = 1000) %>%
    plot_prediction(alpha = 0.65) +
    theme(legend.position = "bottom")
  output_5<-Cross_Validation_lstm(DJ_ordered_new_df,Feature = Features[i],epochs = 1000) %>%
    plot_prediction(alpha = 0.65) +
    theme(legend.position = "bottom")
  output_5<-Cross_Validation_lstm(DM_ordered_new_df,Feature = Features[i],epochs = 1000) %>%
    plot_prediction(alpha = 0.65) +
    theme(legend.position = "bottom")
  pt<-ggarrange(output_1 +   rremove("xlab")+rremove("ylab"), 
                output_2 +   rremove("xlab")+rremove("ylab"), 
                output_3 +   rremove("xlab")+rremove("ylab"), 
                output_4 +   rremove("xlab")+rremove("ylab"), 
                output_5 +   rremove("xlab")+rremove("ylab"),
                output_6 +   rremove("xlab")+rremove("ylab"),
                labels = c("DA", "DB", "DE","DF","DJ","DM"),
                ncol = 6, nrow = 1,common.legend = TRUE,legend="bottom")
  pt<-annotate_figure(pt,top = text_grob(paste0(Features[i]), 
                                         color = "black", face = "bold", size = 14))
  pred_genes[[i]]<-pt
}

