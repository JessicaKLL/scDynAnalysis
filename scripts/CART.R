library(scDynAnalysis)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(pheatmap)
library(ggformula)

DATA<-readRDS("./seurat.rds")
load("./Splited_Patient_data.RData")

#Random Forest

RF_P1<-list()
for (i in 1:4) {
  RF_P1[[i]]<-RandomForest(as.factor(cell_type)~.,data=P1[[i]][,1:1387],strata = "cell_type",ntree=200,split=T)
}
RF_DB<-list()
for (i in 1:4) {
  RF_DB[[i]]<-RandomForest(as.factor(cell_type)~.,data=DB[[i]][,1:1387],strata = "cell_type",ntree=200,split=T)
}
RF_DE<-list()
for (i in 1:4) {
  RF_DE[[i]]<-RandomForest(as.factor(cell_type)~.,data=DE[[i]][,1:1387],strata = "cell_type",ntree=200,split=T)
}
RF_DF<-list()
for (i in 1:4) {
  RF_DF[[i]]<-RandomForest(as.factor(cell_type)~.,data=DF[[i]][,1:1387],strata = "cell_type",ntree=200,split=T)
}
RF_DJ<-list()
for (i in 1:4) {
  RF_DJ[[i]]<-RandomForest(as.factor(cell_type)~.,data=DJ[[i]][,1:1387],strata = "cell_type",ntree=200,split=T)
}

#Quality check

P1_quality<-list()
for (i in 1:length(RF_P1)) {
  P1_quality[[i]]<-model_quality(RF_P1[[i]]$RF_model,RF_P1[[i]]$Test_Set,RF_P1[[i]]$Test_Set$cell_type)
}
names(P1_quality)<-names(P1)
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

P1_impFeat<-list()
for (i in 1:length(P1)) {
  P1_impFeat[[i]]<-important_feat(RF_P1[[i]]$RF_model,data=P1[[i]],n=30, main=paste0(names(P1)[i],": Important Features"))
}
names(P1_impFeat)<-names(P1)
P1_selected_feat<-c()
for (x in 1:length(P1_impFeat)) {
  for (y in 1:length(P1_impFeat[[x]]$Selected)) {
    if (P1_impFeat[[x]]$Selected[y] %in% P1_selected_feat){
      next
    }
    else{
      P1_selected_feat<-append(P1_selected_feat,P1_impFeat[[x]]$Selected[y])
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

gamma_delta_genes<-c("GZMM","GZMK","GZMH","GZMB","GZMA","TYROBP","TRDV1","TRGV8","FCGR3A","KLRF1","NKG7","GNLY","PRF1","CCL4","ZFP36")
exhaustion_markers<-c("PDCD1","LAG3","HAVCR2","TIGIT","CTLA4","ENTPD1","CD244","EOMES","BATF")

Selected_feat<-c(P1_selected_feat,P2_selected_feat,P3_selected_feat,P4_selected_feat,P5_selected_feat,gamma_delta_genes,exhaustion_markers)
Selected_feat<-unique(Selected_feat)

# First filtering

for (i in 1:length(P1)) {
  cell_type<-P1[[i]]$cell_type
  P1[[i]]<-P1[[i]][,which(colnames(P1[[i]]) %in% Selected_feat)]
  P1[[i]]$cell_type<-cell_type
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

# Technical variance signicance

#Patient 1
suppressWarnings(P1_Within_var_IP_Peak<-Within_variance(P1$IP,P1$Peak,n=100,perc=0.3))

#Patient 2
P2$IP$cell_type<-as.factor(P2$IP$cell_type)
P2$Peak$cell_type<-as.factor(P2$Peak$cell_type)
levels(P2$IP$cell_type)<-levels(P2$Peak$cell_type)
suppressWarnings(P2_Within_var_IP_Peak<-Within_variance(P2$IP,P2$Peak,n=100,perc=0.3))

#Patient 3
suppressWarnings(P3_Within_var_IP_Peak<-Within_variance(P3$IP,P3$Peak,n=100,perc=0.3))

#Patient 4
suppressWarnings(P4_Within_var_IP_Peak<-Within_variance(P4$IP,P4$Peak,n=100,perc=0.3))

#Patient 5
suppressWarnings(P5_Within_var_IP_Peak<-Within_variance(P5$IP,P5$Peak,n=100,perc=0.3))

p_values<-data.frame(P1_IP_Peak=mean(P1_Within_var_IP_Peak$p.value),
                     P2_IP_Peak=mean(P2_Within_var_IP_Peak$p.value),
                     P3_IP_Peak=mean(P3_Within_var_IP_Peak$p.value),
                     P4_IP_Peak=mean(P4_Within_var_IP_Peak$p.value),
                     P5_IP_Peak=mean(P5_Within_var_IP_Peak$p.value))

# Altered features

P1_IP<-P1$IP
P1_Peal<-P1$Peak
P1_IP$time_point<-"IP"
P1_Peak$time_point<-"Peak"
P1_complete<-rbind(P1_IP,P1_Peak)
P1_Between_var<-Between_variance(P1_complete,Selected_feat,P1_complete$time_point)
P1_pVal_alter<-select_alter(P1_Between_var,Selected_feat,explicit = T)
P2_IP<-P2$IP
P2_Peal<-P2$Peak
P2_IP$time_point<-"IP"
P2_Peak$time_point<-"Peak"
P2_complete<-rbind(P2_IP,P2_Peak)
P2_Between_var<-Between_variance(P2_complete,Selected_feat,P2_complete$time_point)
P2_pVal_alter<-select_alter(P2_Between_var,Selected_feat,explicit = T)
P3_IP<-P3$IP
P3_Peal<-P3$Peak
P3_IP$time_point<-"IP"
P3_Peak$time_point<-"Peak"
P3_complete<-rbind(P3_IP,P3_Peak)
P3_Between_var<-Between_variance(P3_complete,Selected_feat,P3_complete$time_point)
P3_pVal_alter<-select_alter(P3_Between_var,Selected_feat,explicit = T)
P4_IP<-P4$IP
P4_Peal<-P4$Peak
P4_IP$time_point<-"IP"
P4_Peak$time_point<-"Peak"
P4_complete<-rbind(P4_IP,P4_Peak)
P4_Between_var<-Between_variance(P4_complete,Selected_feat,P4_complete$time_point)
P4_pVal_alter<-select_alter(P4_Between_var,Selected_feat,explicit = T)
P5_IP<-P5$IP
P5_Peal<-P5$Peak
P5_IP$time_point<-"IP"
P5_Peak$time_point<-"Peak"
P5_complete<-rbind(P5_IP,P5_Peak)
P5_Between_var<-Between_variance(P5_complete,Selected_feat,P5_complete$time_point)
P5_pVal_alter<-select_alter(P5_Between_var,Selected_feat,explicit = T)
pVal_alter<-c(P1_pVal_alter,P2_pVal_alter,P3_pVal_alter,P4_pVal_alter,P5_pVal_alter)
pVal_alter<-unique(pVal_alter)

# Second filtering

for (i in 1:length(P1)) {
  cell_type<-P1[[i]]$cell_type
  P1[[i]]<-P1[[i]][,which(colnames(P1[[i]]) %in% pVal_alter)]
  P1[[i]]$cell_type<-cell_type
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

## feature plots
Features<-pVal_alter
P<-FeaturePlot(DATA, features = Features,cols = c("lightgrey","red2"),split.by = "timepoint")

## boxplots

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

P1_IP_cell_type<-P1$IP$cell_type
P1_peak_cell_type<-P1$Peak$cell_type
P2_IP_cell_type<-P2$IP$cell_type
P2_peak_cell_type<-P2$Peak$cell_type
P3_IP_cell_type<-P3$IP$cell_type
P3_peak_cell_type<-P3$Peak$cell_type
P4_IP_cell_type<-P4$IP$cell_type
P4_peak_cell_type<-P4$Peak$cell_type
P5_IP_cell_type<-P5$IP$cell_type
P5_peak_cell_type<-P5$Peak$cell_type

P1$IP$cell_type<-NULL
P1$Peak$cell_type<-NULL
P2$IP$cell_type<-NULL
P2$Peak$cell_type<-NULL
P3$IP$cell_type<-NULL
P3$Peak$cell_type<-NULL
P4$IP$cell_type<-NULL
P4$Peak$cell_type<-NULL
P5$IP$cell_type<-NULL
P5$Peak$cell_type<-NULL

P1$IP$CAR<-NULL
P1$Peak$CAR<-NULL
P2$IP$CAR<-NULL
P2$Peak$CAR<-NULL
P3$IP$CAR<-NULL
P3$Peak$CAR<-NULL
P4$IP$CAR<-NULL
P4$Peak$CAR<-NULL
P5$IP$CAR<-NULL
P5$Peak$CAR<-NULL

P1$IP$time_point<-NULL
P1$Peak$time_point<-NULL
P2$IP$time_point<-NULL
P2$Peak$time_point<-NULL
P3$IP$time_point<-NULL
P3$Peak$time_point<-NULL
P4$IP$time_point<-NULL
P4$Peak$time_point<-NULL
P5$IP$time_point<-NULL
P5$Peak$time_point<-NULL

P1_hclust<-divHier(P1)
P2_hclust<-divHier(P2)
P3_hclust<-divHier(P3)
P4_hclust<-divHier(P4)
P5_hclust<-divHier(P5)

P1_new_df<-gen_new_data(P1_hclust,method = "sum")
P2_new_df<-gen_new_data(P2_hclust,method = "sum")
P3_new_df<-gen_new_data(P3_hclust,method = "sum")
P4_new_df<-gen_new_data(P4_hclust,method = "sum")
P5_new_df<-gen_new_data(P5_hclust,method = "sum")

P1_corr<-find_transition(P1_new_df,cluster = "clusters",feat_num = feat_num)
P2_corr<-find_transition(P2_new_df,cluster = "clusters",feat_num = feat_num)
P3_corr<-find_transition(P3_new_df,cluster = "clusters",feat_num = feat_num)
P4_corr<-find_transition(P4_new_df,cluster = "clusters",feat_num = feat_num)
P5_corr<-find_transition(P5_new_df,cluster = "clusters",feat_num = feat_num)

P1_IP_order<-cell_state_order(P1_new_df[which(P1_new_df$time_point=="IP"),],root = P1_corr$Identical$x)
P1_Peak_order<-cell_state_order(P1_new_df[which(P1_new_df$time_point=="Peak"),],root = P1_corr$Identical$y)
P1_order<-c(rev(P1_IP_order),P1_Peak_order)
P1_ordered_new_df<-P1_new_df[match(P1_order, P1_new_df$clusters),]

P2_IP_order<-cell_state_order(P2_new_df[which(P2_new_df$time_point=="IP"),],root = P2_corr$Identical$x)
P2_Peak_order<-cell_state_order(P2_new_df[which(P2_new_df$time_point=="Peak"),],root = P2_corr$Identical$y)
P2_order<-c(rev(P2_IP_order),P2_Peak_order)
P2_ordered_new_df<-P2_new_df[match(P2_order, P2_new_df$clusters),]

P3_IP_order<-cell_state_order(P3_new_df[which(P3_new_df$time_point=="IP"),],root = P3_corr$Identical$x)
P3_Peak_order<-cell_state_order(P3_new_df[which(P3_new_df$time_point=="Peak"),],root = P3_corr$Identical$y)
P3_order<-c(rev(P3_IP_order),P3_Peak_order)
P3_ordered_new_df<-P3_new_df[match(P3_order, P3_new_df$clusters),]

P4_IP_order<-cell_state_order(P4_new_df[which(P4_new_df$time_point=="IP"),],root = P4_corr$Identical$x)
P4_Peak_order<-cell_state_order(P4_new_df[which(P4_new_df$time_point=="Peak"),],root = P4_corr$Identical$y)
P4_order<-c(rev(P4_IP_order),P4_Peak_order)
P4_ordered_new_df<-P4_new_df[match(P4_order, P4_new_df$clusters),]

P5_IP_order<-cell_state_order(P5_new_df[which(P5_new_df$time_point=="IP"),],root = P5_corr$Identical$x)
P5_Peak_order<-cell_state_order(P5_new_df[which(P5_new_df$time_point=="Peak"),],root = P5_corr$Identical$y)
P5_order<-c(rev(P5_IP_order),P5_Peak_order)
P5_ordered_new_df<-P5_new_df[match(P5_order, P5_new_df$clusters),]

# Dynamics

## spaghetti
gd_donors_new_df<-list(Patient1=P1_ordered_new_df,Patient2=P2_ordered_new_df,Patient3=p3_ordered_new_df,
                            Patient4=P4_ordered_new_df,Patient5=P5_ordered_new_df)
Features<-pVal_alter
Gene_dyn_Xpatient<-list()
for (i in 1:length(Features)) {
  p1<-as.data.frame(gd_patients_new_df$Patient1[,Features[1]])
  p1$order<-gd_patients_new_df$Patient1$clusters
  p1$Patient<-"Patient1"
  names(p1)<-c("Gene","Meta_cells","Patient")
  p2<-as.data.frame(gd_patients_new_df$Patient2[,Features[i]])
  p2$order<-gd_patients_new_df$Patient1$clusters
  p2$Patient<-"Patient2"
  names(p2)<-c("Gene","Meta_cells","Patient")
  p3<-as.data.frame(gd_patients_new_df$Patient3[,Features[i]])
  p3$order<-gd_patients_new_df$Patient1$clusters
  p3$Patient<-"Patient3"
  names(p3)<-c("Gene","Meta_cells","Patient")
  p4<-as.data.frame(gd_patients_new_df$Patient4[,Features[i]])
  p4$order<-gd_patients_new_df$Patient1$clusters
  p4$Patient<-"Patient4"
  names(p4)<-c("Gene","Meta_cells","Patient")
  p5<-as.data.frame(gd_patients_new_df$Patient5[,Features[i]])
  p5$order<-gd_patients_new_df$Patient1$clusters
  p5$Patient<-"Patient5"
  names(p5)<-c("Gene","Meta_cells","Patient")
  df<-rbind(p1,p2)
  df<-rbind(df,p3)
  df<-rbind(df,p4)
  df<-rbind(df,p5)
  df$Meta_cells <- factor(df$Meta_cells, levels=unique(df$Meta_cells))
  df$Patient<-as.factor(df$Patient)
  p<-plot_dyn(df)+
    geom_vline(xintercept = c(200), linetype=2, size = 0.5)+
    geom_text(aes(x=100, label="\nIP", y=10), colour="black") +
    geom_text(aes(x=300, label="\nPeak", y=10), colour="black")+
    scale_color_manual(values = c(Patient1="brown1",Patient2="forestgreen",Patient3="gold",
                                  Patient4="hotpink",Patient5="yellowgreen"))+
    theme(plot.title = element_text(hjust = 0.5))
  Gene_dyn_Xpatient[[i]]<-p
} 

names(Gene_dyn_Xpatient)<-Features

## heatmaps

heatmaps<-list()
heatmaps[[1]]<-Expr_Heatmap(P1_ordered_new_df,Features,main = "P1")
heatmaps[[2]]<-Expr_Heatmap(P2_ordered_new_df,Features,main = "P2")
heatmaps[[3]]<-Expr_Heatmap(P3_ordered_new_df,Features,main = "P3")
heatmaps[[4]]<-Expr_Heatmap(P4_ordered_new_df,Features,main = "P4")
heatmaps[[5]]<-Expr_Heatmap(P5_ordered_new_df,Features,main = "P5")

# Forecasting

Features<-pVal_alter[-which(pVal_alter=="antigen")]

pred_genes<-list()
for (i in 1:length(Features)) {
  output_1<-Cross_Validation_lstm(P1_ordered_new_df,Feature = Features[i],epochs = 1000) %>%
    plot_prediction(alpha = 0.65) +
    theme(legend.position = "bottom")
  output_2<-Cross_Validation_lstm(P2_ordered_new_df,Feature = Features[i],epochs = 1000) %>%
    plot_prediction(alpha = 0.65) +
    theme(legend.position = "bottom")
  output_3<-Cross_Validation_lstm(P3_ordered_new_df,Feature = Features[i],epochs = 1000) %>%
    plot_prediction(alpha = 0.65) +
    theme(legend.position = "bottom")
  output_4<-Cross_Validation_lstm(P4_ordered_new_df,Feature = Features[i],epochs = 1000) %>%
    plot_prediction(alpha = 0.65) +
    theme(legend.position = "bottom")
  output_5<-Cross_Validation_lstm(P5_ordered_new_df,Feature = Features[i],epochs = 1000) %>%
    plot_prediction(alpha = 0.65) +
    theme(legend.position = "bottom")
  pt<-ggarrange(output_1 +   rremove("xlab")+rremove("ylab"), 
                output_2 +   rremove("xlab")+rremove("ylab"), 
                output_3 +   rremove("xlab")+rremove("ylab"), 
                output_4 +   rremove("xlab")+rremove("ylab"), 
                output_5 +   rremove("xlab")+rremove("ylab"), 
                labels = c("P1", "P2", "P3","P4","P5"),
                ncol = 5, nrow = 1,common.legend = TRUE,legend="bottom")
  pt<-annotate_figure(pt,top = text_grob(paste0(Features[i]), 
                                         color = "black", face = "bold", size = 14))
  pred_genes[[i]]<-pt
}

