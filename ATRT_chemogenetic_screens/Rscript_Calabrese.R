library(MAGeCKFlute)
library(ggplot2)
library(tidyverse)
library(dplyr)

#generate custom function for Scatter plot with SD=2
ScatterViewCustom<-function(data, x = "x", y = "y", label = 0,
                            model = c("none", "ninesquare", "volcano", "rank")[1],
                            x_cut = NULL, y_cut = NULL, slope = 1, intercept = NULL,
                            auto_cut = FALSE, auto_cut_x = auto_cut,
                            auto_cut_y = auto_cut, auto_cut_diag = auto_cut,
                            groups = NULL, group_col = NULL, groupnames = NULL,
                            label.top = TRUE, top = 0, toplabels = NULL,
                            display_cut = FALSE, color = NULL, shape = 16, size = 1, alpha = 0.6,
                            main = NULL, xlab = x, ylab = y, legend.position = "none", ...){
  requireNamespace("ggplot2", quietly=TRUE) || stop("need ggplot package")
  requireNamespace("ggrepel", quietly=TRUE) || stop("need ggrepel package")
  data = as.data.frame(data, stringsAsFactors = FALSE)
  data = data[!(is.na(data[,x])|is.na(data[,y])), ]
  ## Add label column in the data frame.
  if(label==0) data$Label = rownames(data)
  else data$Label = as.character(data[, label])
  
  if(!is.null(groupnames)) legend.position = "right"
  
  ## Compute the cutoff used for each dimension.
  model = tolower(model)
  if(model == "ninesquare"){
    if(length(x_cut)==0)
      x_cut = c(-CutoffCalling(data[,x], 2), CutoffCalling(data[,x], 2))
    if(length(y_cut)==0)
      y_cut = c(-CutoffCalling(data[,y], 2), CutoffCalling(data[,y], 2))
    if(length(intercept)==0)
      intercept = c(-CutoffCalling(data[,y]-data[,x], 2), CutoffCalling(data[,y]-data[,x], 2))
  }
  if(model == "volcano"){
    if(length(x_cut)==0)
      x_cut = c(-CutoffCalling(data[,x], 2), CutoffCalling(data[,x], 2))
    if(length(y_cut)==0) y_cut = -log10(0.05)
  }
  if(model == "rank"){
    if(length(x_cut)==0)
      x_cut = c(-CutoffCalling(data[,x], 2), CutoffCalling(data[,x], 2))
  }
  if(model == "none"){
    if(auto_cut_x)
      x_cut = c(-CutoffCalling(data[,x], 2), CutoffCalling(data[,x], 2))
    if(auto_cut_y)
      y_cut = c(-CutoffCalling(data[,y], 2), CutoffCalling(data[,y], 2))
    if(auto_cut_diag)
      intercept = c(-CutoffCalling(data[,y]-data[,x], 2), CutoffCalling(data[,y]-data[,x], 2))
  }
  ## Decide the colored groups
  avail_groups = c("topleft", "topright", "bottomleft", "bottomright",
                   "midleft", "topcenter", "midright", "bottomcenter", "midcenter",
                   "top", "mid", "bottom", "left", "center", "right", "none")
  ## Select the colors
  mycolour = c("#1f78b4", "#fb8072", "#33a02c", "#ff7f00",
               "#571782", "#e41a1c", "#178261", "#377eb8", "#ffed6f",
               "#e78ac3", "#fdb462", "#8da0cb", "#66c2a5", "#fccde5", "#fc8d62", "#d9d9d9")
  names(mycolour) = avail_groups
  
  if(model == "ninesquare") groups = c("midleft", "topcenter", "midright", "bottomcenter")
  if(model == "volcano") groups = c("topleft", "topright")
  if(model == "rank") groups = c("left", "right")
  groups = intersect(groups, avail_groups)
  
  ## Annotate the groups in the data frame
  if(length(x_cut)>0){
    idx1 = data[,x] < min(x_cut)
    idx2 = data[,x] > max(x_cut)
  }else{
    idx1 = NA
    idx2 = NA
  }
  if(length(y_cut)>0){
    idx3 = data[,y] < min(y_cut)
    idx4 = data[,y] > max(y_cut)
  }else{
    idx3 = NA
    idx4 = NA
  }
  if(length(intercept)>0){
    idx5 = data[,y]<slope*data[,x]+min(intercept)
    idx6 = data[,y]>slope*data[,x]+max(intercept)
  }else{
    idx5 = NA; idx6 = NA
  }
  data$group="none"
  for(gr in groups){
    if(gr=="topleft") idx = cbind(idx1, idx4, idx6)
    if(gr=="topcenter") idx = cbind(!idx1, !idx2, idx4, idx6)
    if(gr=="topright") idx = cbind(idx2, idx4, idx6)
    if(gr=="midleft") idx = cbind(idx1, idx6 , !idx3, !idx4)
    if(gr=="midcenter") idx = cbind(!idx1, !idx2, !idx3, !idx4, !idx5, !idx6)
    if(gr=="midright") idx = cbind(idx2, !idx3, !idx4, idx5)
    if(gr=="bottomleft") idx = cbind(idx1, idx3, idx5)
    if(gr=="bottomcenter") idx = cbind(!idx1, !idx2, idx3, idx5)
    if(gr=="bottomright") idx = cbind(idx2, idx3, idx5)
    if(gr=="top"){
      if(length(y_cut)>0 & length(intercept)>0)
        idx = idx4 & idx6
      else if(length(y_cut)>0)
        idx = idx4
      else idx = idx6
    }
    if(gr=="mid") idx = (!idx3) & (!idx4)
    if(gr=="bottom"){
      if(length(y_cut)>0 & length(intercept)>0)
        idx = idx3 & idx5
      else if(length(y_cut)>0)
        idx = idx3
      else idx = idx5
    }
    if(gr=="left"){
      if(length(x_cut)>0 & length(intercept)>0)
        if(slope>0) idx = idx1 & idx6 else idx = idx1 & idx5
        else if(length(x_cut)>0)
          idx = idx1
        else
          if(slope>0) idx = idx6 else idx = idx5
    }
    if(gr=="center") idx = (!idx1) & (!idx2)
    if(gr=="right"){
      if(length(x_cut)>0 & length(intercept)>0)
        if(slope>0) idx = idx2 & idx5 else idx = idx2 & idx6
        else if(length(x_cut)>0)
          idx = idx2
        else
          if(slope>0) idx = idx5 else idx = idx6
    }
    ## Assign groups
    if(is.null(ncol(idx))){
      if(sum(!is.na(idx))>0) data$group[idx] = gr
      else warning("No cutpoint for group:", gr)
    }else{
      idx = idx[, !is.na(idx[1,])]
      if(is.null(ncol(idx)))
        warning("No cutpoint for group:", gr)
      else if(ncol(idx)<4 & gr=="midcenter")
        warning("No cutpoint for group:", gr)
      else
        data$group[rowSums(idx)==ncol(idx)] = gr
    }
  }
  data$group=factor(data$group, levels = unique(c(groups, "none")))
  ## Group names
  if(length(groupnames)!=length(groups)) groupnames = groups
  if(length(groups)>0) names(groupnames) = groups
  if(length(group_col)==length(groups)) mycolour[groups] = group_col
  if(length(groups)==0) mycolour["none"] = "#FF6F61"
  
  ## Label top gene names ##
  data$rank = top + 1
  for(g in groups){
    idx1 = data$group==g
    x_symb = 0; y_symb = 0;
    if(g=="topleft"){ x_symb = 1; y_symb = -1 }
    if(g=="topcenter"){ x_symb = 0; y_symb = -1 }
    if(g=="topright"){ x_symb = -1; y_symb = -1 }
    if(g=="midleft"){ x_symb = 1; y_symb = 0 }
    if(g=="midright"){ x_symb = -1; y_symb = 0 }
    if(g=="bottomleft"){ x_symb = 1; y_symb = 1 }
    if(g=="bottomcenter"){ x_symb = 0; y_symb = 1 }
    if(g=="bottomright"){ x_symb = -1; y_symb = 1 }
    if(g=="top"){ x_symb = 0; y_symb = -1 }
    if(g=="bottom"){ x_symb = 0; y_symb = 1 }
    if(g=="left"){ x_symb = 1; y_symb = 0 }
    if(g=="right"){ x_symb = -1; y_symb = 0 }
    tmp = data[,c(x,y)]
    tmp[,x] = (tmp[,x]-min(tmp[,x])) / (max(tmp[,x])-min(tmp[,x]))
    tmp[,y] = (tmp[,y]-min(tmp[,y])) / (max(tmp[,y])-min(tmp[,y]))
    data$rank[idx1] = rank((x_symb*tmp[,x]+y_symb*tmp[,y])[idx1])
  }
  data$rank[data$rank==0] = Inf
  if(mode(toplabels)=="list"){
    data$Label[data$rank>top & !(data$Label %in% unlist(toplabels))] = ""
    data$group = data$Label;
    if(length(toplabels)>0){
      tmp = stack(toplabels)
      tmp = tmp[!duplicated(tmp[,1]), ]
      rownames(tmp) = tmp[,1]
      data$group[data$group%in%tmp[,1]] = as.character(tmp[data$group[data$group%in%tmp[,1]], 2])
      data$group[!(data$group%in%tmp[,2]) & data$group!=""] = "Top hits"
    }
  }else{
    data$Label[data$rank>top & !(data$Label %in% toplabels)] = ""
  }
  
  ## Color issue
  if(is.null(color)){
    color = "group"
  }else if(length(color)==1){
    if(!color%in%colnames(data)){
      data$color = color
      color = "color"
    }
  }else{
    data$color = color[1]
    color = "color"
    warning("Only the first color is took.")
  }
  
  ## Plot the scatter figure ##
  gg = data
  
  ## Plot the figure
  gg = gg[order(gg[,color]), ]
  p = ggplot(gg, aes_string(x, y, label="Label", color = color))
  if(all(c(shape,size)%in%colnames(gg)))
    p = p + geom_point(aes_string(shape = shape, size = size), alpha = alpha)
  else if(shape%in%colnames(gg))
    p = p + geom_point(aes_string(shape = shape), size = size, alpha = alpha)
  else if(size%in%colnames(gg))
    p = p + geom_point(aes_string(size = size), shape = shape, alpha = alpha)
  else
    p = p + geom_point(size = size, shape = shape, alpha = alpha)
  
  ## Customize colors
  if(color=="group"){
    if(mode(toplabels)!="list")
      p = p + scale_color_manual(values = mycolour, labels = groupnames)
    else
      p = p + scale_color_manual(values = c("#d9d9d9", "#fb8072", "#80b1d3", "#fdb462",
                                            "#bc80bd", "#b3de69", "#bebada", "#8dd3c7",
                                            "#ffffb3", "#fccde5", "#ccebc5", "#ffed6f"))
  }else{
    if(mode(gg[,color])=="numeric")
      p = p + scale_color_gradient2(low = "#377eb8", high = "#e41a1c", midpoint = 0)
    else if(!"try-error"%in%class(try(col2rgb(gg[1,color]),silent=TRUE))){
      mycolour = unique(gg[,color]); names(mycolour) = mycolour
      p = p + scale_color_manual(values = mycolour)
    }else{
      p = p + scale_color_brewer(type = "div")
    }
  }
  
  if(label.top)
    p = p + ggrepel::geom_text_repel(...)
  if(display_cut){
    if(length(x_cut)>0)
      p = p + geom_vline(xintercept = x_cut,linetype = "dotted")
    if(length(y_cut)>0)
      p = p + geom_hline(yintercept = y_cut,linetype = "dotted")
    if(length(intercept)>0)
      p = p + geom_abline(slope=slope, intercept=intercept, linetype = "dotted")
  }
  p = p + labs(x=xlab, y = ylab, title = main, color = NULL)
  p = p + theme_bw(base_size = 12)
  p = p + theme(legend.position = legend.position)
  
  return(p)
}

###Analysis and visualization of Brunello screens as drug common changes in BT16 and CHLA06
###work with MLE output from design drug_common vs plasmid and DMSO vs plasmid per line for Scatter, drug_common vs DMSO for Rank plot
#get data
getwd()
setwd("\\\\fileserv1/groups/Neuro-Oncology/Clinical Neuro-Oncology/Daniel/B22-ATRT screens/Modifier_screens/Calabrese/MAGECK Flute/Final_Analysis/data")

all.screens.vs.plasmid = read.delim("All_Screens_vs_plasmid.gene_summary.txt")
all.screens.vs.plasmid = all.screens.vs.plasmid[-129,]
data.BT16.DMSOvsplasmid = all.screens.vs.plasmid[,c(1,3)]
data.CHLA06.DMSOvsplasmid = all.screens.vs.plasmid[,-(3:20)]
data.CHLA06.DMSOvsplasmid = data.CHLA06.DMSOvsplasmid[,1:8]
data.CHLA06.DMSOvsplasmid = data.CHLA06.DMSOvsplasmid[,c(1,3)]

data.BT16.drugcommonvsplasmid = ReadBeta("BT16_drug_common_vs_plasmid.gene_summary.txt")
data.BT16.scatter = cbind(data.BT16.drugcommonvsplasmid, data.BT16.DMSOvsplasmid$BT16_DMSO.beta)
colnames(data.BT16.scatter)[colnames(data.BT16.scatter) == 'BT16_drug_common_vs_plasmid'] <- 'CDKi'
colnames(data.BT16.scatter)[colnames(data.BT16.scatter) == 'data.BT16.DMSOvsplasmid$BT16_DMSO.beta'] <- 'DMSO'
data.BT16.DrugvsDMSO =ReadBeta("BT16_drug_common_vs_DMSO.gene_summary.txt")
data.BT16.DrugvsDMSO = data.BT16.DrugvsDMSO[-129,]
colnames(data.BT16.DrugvsDMSO)[colnames(data.BT16.DrugvsDMSO) == 'BT16_drug_common'] <- 'CDKi_vs_DMSO'

data.CHLA06.drugcommonvsplasmid = ReadBeta("CHLA06_drug_common_vs_plasmid.gene_summary.txt")
data.CHLA06.scatter = cbind(data.CHLA06.drugcommonvsplasmid, data.CHLA06.DMSOvsplasmid$CHLA06_DMSO.beta)
colnames(data.CHLA06.scatter)[colnames(data.CHLA06.scatter) == 'CHLA06_drug_common_vs_plasmid'] <- 'CDKi'
colnames(data.CHLA06.scatter)[colnames(data.CHLA06.scatter) == 'data.CHLA06.DMSOvsplasmid$CHLA06_DMSO.beta'] <- 'DMSO'
data.CHLA06.DrugvsDMSO =ReadBeta("CHLA06_drug_common_vs_DMSO.gene_summary.txt")
colnames(data.CHLA06.DrugvsDMSO)[colnames(data.CHLA06.DrugvsDMSO) == 'CHLA06_drug_common'] <- 'CDKi_vs_DMSO'


#make Scatter and rank plot for BT16 cells, in rank plot show selected genes with FDR < 0.1 and those also in scatter
setwd("\\\\fileserv1/groups/Neuro-Oncology/Clinical Neuro-Oncology/Daniel/B22-ATRT screens/Modifier_screens/Calabrese/MAGECK Flute/Final_Analysis/results")
options(ggrepel.max.overlaps = Inf)
genelistBT16.CDKivsDMSO = data.BT16.DrugvsDMSO$CDKi_vs_DMSO
names(genelist.BT16.CDKivsDMSO) =data.BT16.DrugvsDMSO$Gene
Rank.BT16.CDKivsDMSO = RankView(genelist.BT16.CDKivsDMSO,
                                      top = 5, bottom = 5,
                                      main = "Rank BT16 CDKivsDMSO",
                                      width = 2, height = 4)
print(Rank.BT16.CDKivsDMSO + theme_classic()+theme(legend.position="none"), vp=grid::viewport(gp=grid::gpar(cex=2)))

Rank.BT16.CDKivsDMSO.clean = RankView(genelist.BT16.CDKivsDMSO,
                                      top = 0, bottom = 0,
                                      main = "Rank BT16 CDKivsDMSO",
                                      width = 2, height = 4)
print(Rank.BT16.CDKivsDMSO.clean + theme_classic()+theme(legend.position="none"), vp=grid::viewport(gp=grid::gpar(cex=2)))

Scatter.BT16.CDKi =ScatterViewCustom(data.BT16.scatter, x= "DMSO", y= "CDKi", 
                                     label = "Gene",size =2,groups = c("top", "bottom"), top=0,
                                     toplabels = c( "PLAUR", "DSG4", "RAET1L", "POTEG"),
                                     auto_cut_diag  = 2, group_col = c("#e41a1c", "#377eb8"),
                                     alpha = 0.5, display_cut = TRUE)
print(Scatter.BT16.CDKi + theme_classic() + theme(legend.position="none"))

Scatter.BT16.CDKi.9sq =ScatterViewCustom(data.BT16.scatter, x= "DMSO", y= "CDKi", 
                                         label = "Gene",size =2,groups = c("top", "bottom"),
                                         model = c("ninesquare"),
                                         auto_cut_diag  = 2, group_col = c("#e41a1c", "#377eb8"),
                                         alpha = 0.5, display_cut = TRUE)
print(Scatter.BT16.CDKi.9sq + theme_classic() + theme(legend.position="none"))

Scatter.BT16.CDKi.clean =ScatterViewCustom(data.BT16.scatter, x= "DMSO", y= "CDKi", 
                                           label = "Gene",size =2,groups = c("top", "bottom"),
                                           auto_cut_diag  = 2, group_col = c("#e41a1c", "#377eb8"),
                                           alpha = 0.5, display_cut = TRUE)
print(Scatter.BT16.CDKi.clean + theme_classic() + theme(legend.position="none"))


#make Scatter and rank plot for CHLA06 cells, in rank plot show selected genes with FDR < 0.1 and those also in scatter
options(ggrepel.max.overlaps = Inf)

genelist.CHLA06.CDKivsDMSO = data.CHLA06.DrugvsDMSO$CDKi_vs_DMSO
names(genelist.CHLA06.CDKivsDMSO) =data.CHLA06.DrugvsDMSO$Gene
Rank.CHLA06.CDKivsDMSO = RankView(genelist.CHLA06.CDKivsDMSO,
                                  top = 5, bottom = 5,
                                  genelist = c("AMBRA1", "FBXW7", "UBA1"),
                                  main = "Rank CHLA06 CDKivsDMSO",
                                  width = 2, height = 4)
print(Rank.CHLA06.CDKivsDMSO + theme_classic()+theme(legend.position="none"), vp=grid::viewport(gp=grid::gpar(cex=2)))


Rank.CHLA06.CDKivsDMSO.clean = RankView(genelist.CHLA06.CDKivsDMSO,
                                        top = 0, bottom = 0,
                                        main = "Rank CHLA06 CDKivsDMSO",
                                        width = 2, height = 4)
print(Rank.CHLA06.CDKivsDMSO.clean + theme_classic()+theme(legend.position="none"), vp=grid::viewport(gp=grid::gpar(cex=2)))

Scatter.CHLA06.CDKi =ScatterViewCustom(data.CHLA06.scatter, x= "DMSO", y= "CDKi", 
                                       label = "Gene",size =2,groups = c("top", "bottom"), top=0,
                                       toplabels = c("CCNE1", "CCNE2", "KIF2C", "SMC2", "POTEG", "CCNG2"),
                                       auto_cut_diag  = 2, group_col = c("#e41a1c", "#377eb8"),
                                       alpha = 0.5, display_cut = TRUE)
print(Scatter.CHLA06.CDKi + theme_classic() + theme(legend.position="none"))

Scatter.CHLA06.CDKi.9sq =ScatterViewCustom(data.CHLA06.scatter, x= "DMSO", y= "CDKi", 
                                           label = "Gene",size =2,groups = c("top", "bottom"),
                                           model = c("ninesquare"),
                                           auto_cut_diag  = 2, group_col = c("#e41a1c", "#377eb8"),
                                           alpha = 0.5, display_cut = TRUE)
print(Scatter.CHLA06.CDKi.9sq + theme_classic() + theme(legend.position="none"))

Scatter.CHLA06.CDKi.clean =ScatterViewCustom(data.CHLA06.scatter, x= "DMSO", y= "CDKi", 
                                             label = "Gene",size =2,groups = c("top", "bottom"),
                                             auto_cut_diag  = 2, group_col = c("#e41a1c", "#377eb8"),
                                             alpha = 0.5, display_cut = TRUE)
print(Scatter.CHLA06.CDKi.clean + theme_classic() + theme(legend.position="none"))
