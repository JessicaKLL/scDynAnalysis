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
for (i in 1:2) {
  RF_P1[[i]]<-RandomForest(as.factor(cell_type)~.,data=P1[[i]][,1:1387],strata = "cell_type",ntree=200,split=T)
}
RF_P2<-list()
for (i in 1:2) {
  RF_P2[[i]]<-RandomForest(as.factor(cell_type)~.,data=P2[[i]][,1:1387],strata = "cell_type",ntree=200,split=T)
}
RF_P3<-list()
for (i in 1:2) {
  RF_P3[[i]]<-RandomForest(as.factor(cell_type)~.,data=P3[[i]][,1:1387],strata = "cell_type",ntree=200,split=T)
}
RF_P4<-list()
for (i in 1:2) {
  RF_P4[[i]]<-RandomForest(as.factor(cell_type)~.,data=P4[[i]][,1:1387],strata = "cell_type",ntree=200,split=T)
}
RF_P5<-list()
for (i in 1:2) {
  RF_P5[[i]]<-RandomForest(as.factor(cell_type)~.,data=P5[[i]][,1:1387],strata = "cell_type",ntree=200,split=T)
}

#Quality check

P1_quality<-list()
for (i in 1:length(RF_P1)) {
  P1_quality[[i]]<-model_quality(RF_P1[[i]]$RF_model,RF_P1[[i]]$Test_Set,RF_P1[[i]]$Test_Set$cell_type)
}
names(P1_quality)<-names(P1)
P2_quality<-list()
for (i in 1:length(RF_P2)) {
  P2_quality[[i]]<-model_quality(RF_P2[[i]]$RF_model,RF_P2[[i]]$Test_Set,RF_P2[[i]]$Test_Set$cell_type)
}
names(P2_quality)<-names(P2)
P3_quality<-list()
for (i in 1:length(RF_P3)) {
  P3_quality[[i]]<-model_quality(RF_P3[[i]]$RF_model,RF_P3[[i]]$Test_Set,RF_P3[[i]]$Test_Set$cell_type)
}
names(P3_quality)<-names(P3)
P4_quality<-list()
for (i in 1:length(RF_P4)) {
  P4_quality[[i]]<-model_quality(RF_P4[[i]]$RF_model,RF_P4[[i]]$Test_Set,RF_P4[[i]]$Test_Set$cell_type)
}
names(P4_quality)<-names(P4)
P5_quality<-list()
for (i in 1:length(RF_P5)) {
  P5_quality[[i]]<-model_quality(RF_P5[[i]]$RF_model,RF_P5[[i]]$Test_Set,RF_P5[[i]]$Test_Set$cell_type)
}
names(P5_quality)<-names(P5)

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
P2_impFeat<-list()
for (i in 1:length(P2)) {
  P2_impFeat[[i]]<-important_feat(RF_P2[[i]]$RF_model,data=P2[[i]],n=30, main=paste0(names(P2)[i],": Important Features"))
}
names(P2_impFeat)<-names(P2)
P2_selected_feat<-c()
for (x in 1:length(P2_impFeat)) {
  for (y in 1:length(P2_impFeat[[x]]$Selected)) {
    if (P2_impFeat[[x]]$Selected[y] %in% P2_selected_feat){
      next
    }
    else{
      P2_selected_feat<-append(P2_selected_feat,P2_impFeat[[x]]$Selected[y])
    }
  }
}
P3_impFeat<-list()
for (i in 1:length(P3)) {
  P3_impFeat[[i]]<-important_feat(RF_P3[[i]]$RF_model,data=P3[[i]],n=30, main=paste0(names(P3)[i],": Important Features"))
}
names(P3_impFeat)<-names(P3)
P3_selected_feat<-c()
for (x in 1:length(P3_impFeat)) {
  for (y in 1:length(P3_impFeat[[x]]$Selected)) {
    if (P3_impFeat[[x]]$Selected[y] %in% P3_selected_feat){
      next
    }
    else{
      P3_selected_feat<-append(P3_selected_feat,P3_impFeat[[x]]$Selected[y])
    }
  }
}
P4_impFeat<-list()
for (i in 1:length(P4)) {
  P4_impFeat[[i]]<-important_feat(RF_P4[[i]]$RF_model,data=P4[[i]],n=30, main=paste0(names(P4)[i],": Important Features"))
}
names(P4_impFeat)<-names(P4)
P4_selected_feat<-c()
for (x in 1:length(P4_impFeat)) {
  for (y in 1:length(P4_impFeat[[x]]$Selected)) {
    if (P4_impFeat[[x]]$Selected[y] %in% P4_selected_feat){
      next
    }
    else{
      P4_selected_feat<-append(P4_selected_feat,P4_impFeat[[x]]$Selected[y])
    }
  }
}
P5_impFeat<-list()
for (i in 1:length(P5)) {
  P5_impFeat[[i]]<-important_feat(RF_P5[[i]]$RF_model,data=P5[[i]],n=30, main=paste0(names(P5)[i],": Important Features"))
}
names(P5_impFeat)<-names(P5)
P5_selected_feat<-c()
for (x in 1:length(P5_impFeat)) {
  for (y in 1:length(P5_impFeat[[x]]$Selected)) {
    if (P5_impFeat[[x]]$Selected[y] %in% P5_selected_feat){
      next
    }
    else{
      P5_selected_feat<-append(P5_selected_feat,P5_impFeat[[x]]$Selected[y])
    }
  }
}

gamma_P3lta_genes<-c("GZMM","GZMK","GZMH","GZMB","GZMA","TYROBP","TRDV1","TRGV8","FCGR3A","KLRF1","NKG7","GNLY","PRF1","CCL4","ZFP36")
exhaustion_markers<-c("PDCD1","LAG3","HAVCR2","TIGIT","CTLA4","ENTPD1","CD244","EOMES","BATF")

Selected_feat<-c(P1_selected_feat,P2_selected_feat,P3_selected_feat,P4_selected_feat,P5_selected_feat,gamma_P3lta_genes,exhaustion_markers)
Selected_feat<-unique(Selected_feat)

# First filtering

for (i in 1:length(P1)) {
  cell_type<-P1[[i]]$cell_type
  P1[[i]]<-P1[[i]][,which(colnames(P1[[i]]) %in% Selected_feat)]
  P1[[i]]$cell_type<-cell_type
}
for (i in 1:length(P2)) {
  cell_type<-P2[[i]]$cell_type
  P2[[i]]<-P2[[i]][,which(colnames(P2[[i]]) %in% Selected_feat)]
  P2[[i]]$cell_type<-cell_type
}
for (i in 1:length(P3)) {
  cell_type<-P3[[i]]$cell_type
  P3[[i]]<-P3[[i]][,which(colnames(P3[[i]]) %in% Selected_feat)]
  P3[[i]]$cell_type<-cell_type
}
for (i in 1:length(P4)) {
  cell_type<-P4[[i]]$cell_type
  P4[[i]]<-P4[[i]][,which(colnames(P4[[i]]) %in% Selected_feat)]
  P4[[i]]$cell_type<-cell_type
}
for (i in 1:length(P5)) {
  cell_type<-P5[[i]]$cell_type
  P5[[i]]<-P5[[i]][,which(colnames(P5[[i]]) %in% Selected_feat)]
  P5[[i]]$cell_type<-cell_type
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
for (i in 1:length(P2)) {
  cell_type<-P2[[i]]$cell_type
  P2[[i]]<-P2[[i]][,which(colnames(P2[[i]]) %in% pVal_alter)]
  P2[[i]]$cell_type<-cell_type
}
for (i in 1:length(P3)) {
  cell_type<-P3[[i]]$cell_type
  P3[[i]]<-P3[[i]][,which(colnames(P3[[i]]) %in% pVal_alter)]
  P3[[i]]$cell_type<-cell_type
}
for (i in 1:length(P4)) {
  cell_type<-P4[[i]]$cell_type
  P4[[i]]<-P4[[i]][,which(colnames(P4[[i]]) %in% pVal_alter)]
  P4[[i]]$cell_type<-cell_type
}
for (i in 1:length(P5)) {
  cell_type<-P5[[i]]$cell_type
  P5[[i]]<-P5[[i]][,which(colnames(P5[[i]]) %in% pVal_alter)]
  P5[[i]]$cell_type<-cell_type
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

# cell state orP3ring

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

P1_new_P4<-gen_new_data(P1_hclust,method = "sum")
P2_new_P4<-gen_new_data(P2_hclust,method = "sum")
P3_new_P4<-gen_new_data(P3_hclust,method = "sum")
P4_new_P4<-gen_new_data(P4_hclust,method = "sum")
P5_new_P4<-gen_new_data(P5_hclust,method = "sum")

P1_corr<-find_transition(P1_new_P4,cluster = "clusters",feat_num = feat_num)
P2_corr<-find_transition(P2_new_P4,cluster = "clusters",feat_num = feat_num)
P3_corr<-find_transition(P3_new_P4,cluster = "clusters",feat_num = feat_num)
P4_corr<-find_transition(P4_new_P4,cluster = "clusters",feat_num = feat_num)
P5_corr<-find_transition(P5_new_P4,cluster = "clusters",feat_num = feat_num)

P1_IP_orP3r<-cell_state_orP3r(P1_new_P4[which(P1_new_P4$time_point=="IP"),],root = P1_corr$IP3ntical$x)
P1_Peak_orP3r<-cell_state_orP3r(P1_new_P4[which(P1_new_P4$time_point=="Peak"),],root = P1_corr$IP3ntical$y)
P1_orP3r<-c(rev(P1_IP_orP3r),P1_Peak_orP3r)
P1_orP3red_new_P4<-P1_new_P4[match(P1_orP3r, P1_new_P4$clusters),]

P2_IP_orP3r<-cell_state_orP3r(P2_new_P4[which(P2_new_P4$time_point=="IP"),],root = P2_corr$IP3ntical$x)
P2_Peak_orP3r<-cell_state_orP3r(P2_new_P4[which(P2_new_P4$time_point=="Peak"),],root = P2_corr$IP3ntical$y)
P2_orP3r<-c(rev(P2_IP_orP3r),P2_Peak_orP3r)
P2_orP3red_new_P4<-P2_new_P4[match(P2_orP3r, P2_new_P4$clusters),]

P3_IP_orP3r<-cell_state_orP3r(P3_new_P4[which(P3_new_P4$time_point=="IP"),],root = P3_corr$IP3ntical$x)
P3_Peak_orP3r<-cell_state_orP3r(P3_new_P4[which(P3_new_P4$time_point=="Peak"),],root = P3_corr$IP3ntical$y)
P3_orP3r<-c(rev(P3_IP_orP3r),P3_Peak_orP3r)
P3_orP3red_new_P4<-P3_new_P4[match(P3_orP3r, P3_new_P4$clusters),]

P4_IP_orP3r<-cell_state_orP3r(P4_new_P4[which(P4_new_P4$time_point=="IP"),],root = P4_corr$IP3ntical$x)
P4_Peak_orP3r<-cell_state_orP3r(P4_new_P4[which(P4_new_P4$time_point=="Peak"),],root = P4_corr$IP3ntical$y)
P4_orP3r<-c(rev(P4_IP_orP3r),P4_Peak_orP3r)
P4_orP3red_new_P4<-P4_new_P4[match(P4_orP3r, P4_new_P4$clusters),]

P5_IP_orP3r<-cell_state_orP3r(P5_new_P4[which(P5_new_P4$time_point=="IP"),],root = P5_corr$IP3ntical$x)
P5_Peak_orP3r<-cell_state_orP3r(P5_new_P4[which(P5_new_P4$time_point=="Peak"),],root = P5_corr$IP3ntical$y)
P5_orP3r<-c(rev(P5_IP_orP3r),P5_Peak_orP3r)
P5_orP3red_new_P4<-P5_new_P4[match(P5_orP3r, P5_new_P4$clusters),]

# Dynamics

## spaghetti
gd_donors_new_P4<-list(Patient1=P1_orP3red_new_P4,Patient2=P2_orP3red_new_P4,Patient3=p3_orP3red_new_P4,
                            Patient4=P4_orP3red_new_P4,Patient5=P5_orP3red_new_P4)
Features<-pVal_alter
Gene_dyn_Xpatient<-list()
for (i in 1:length(Features)) {
  p1<-as.data.frame(gd_patients_new_P4$Patient1[,Features[1]])
  p1$orP3r<-gd_patients_new_P4$Patient1$clusters
  p1$Patient<-"Patient1"
  names(p1)<-c("Gene","Meta_cells","Patient")
  p2<-as.data.frame(gd_patients_new_P4$Patient2[,Features[i]])
  p2$orP3r<-gd_patients_new_P4$Patient1$clusters
  p2$Patient<-"Patient2"
  names(p2)<-c("Gene","Meta_cells","Patient")
  p3<-as.data.frame(gd_patients_new_P4$Patient3[,Features[i]])
  p3$orP3r<-gd_patients_new_P4$Patient1$clusters
  p3$Patient<-"Patient3"
  names(p3)<-c("Gene","Meta_cells","Patient")
  p4<-as.data.frame(gd_patients_new_P4$Patient4[,Features[i]])
  p4$orP3r<-gd_patients_new_P4$Patient1$clusters
  p4$Patient<-"Patient4"
  names(p4)<-c("Gene","Meta_cells","Patient")
  p5<-as.data.frame(gd_patients_new_P4$Patient5[,Features[i]])
  p5$orP3r<-gd_patients_new_P4$Patient1$clusters
  p5$Patient<-"Patient5"
  names(p5)<-c("Gene","Meta_cells","Patient")
  P4<-rbind(p1,p2)
  P4<-rbind(P4,p3)
  P4<-rbind(P4,p4)
  P4<-rbind(P4,p5)
  P4$Meta_cells <- factor(P4$Meta_cells, levels=unique(P4$Meta_cells))
  P4$Patient<-as.factor(P4$Patient)
  p<-plot_dyn(P4)+
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
heatmaps[[1]]<-Expr_Heatmap(P1_orP3red_new_P4,Features,main = "P1")
heatmaps[[2]]<-Expr_Heatmap(P2_orP3red_new_P4,Features,main = "P2")
heatmaps[[3]]<-Expr_Heatmap(P3_orP3red_new_P4,Features,main = "P3")
heatmaps[[4]]<-Expr_Heatmap(P4_orP3red_new_P4,Features,main = "P4")
heatmaps[[5]]<-Expr_Heatmap(P5_orP3red_new_P4,Features,main = "P5")

# Forecasting

Features<-pVal_alter[-which(pVal_alter=="antigen")]

pred_genes<-list()
for (i in 1:length(Features)) {
  output_1<-Cross_Validation_lstm(P1_orP3red_new_P4,Feature = Features[i],epochs = 1000) %>%
    plot_prediction(alpha = 0.65) +
    theme(legend.position = "bottom")
  output_2<-Cross_Validation_lstm(P2_orP3red_new_P4,Feature = Features[i],epochs = 1000) %>%
    plot_prediction(alpha = 0.65) +
    theme(legend.position = "bottom")
  output_3<-Cross_Validation_lstm(P3_orP3red_new_P4,Feature = Features[i],epochs = 1000) %>%
    plot_prediction(alpha = 0.65) +
    theme(legend.position = "bottom")
  output_4<-Cross_Validation_lstm(P4_orP3red_new_P4,Feature = Features[i],epochs = 1000) %>%
    plot_prediction(alpha = 0.65) +
    theme(legend.position = "bottom")
  output_5<-Cross_Validation_lstm(P5_orP3red_new_P4,Feature = Features[i],epochs = 1000) %>%
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

