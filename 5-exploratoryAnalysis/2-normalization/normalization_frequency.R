################################################################################
## normalization in tcr and bcr
## jean resende 
# big data
################################################################################
# -- normalization in tcr and bcr: frequency column -- #########################
## -- region V -- 
load("../1-freeze_date/freeze_tcrbcr_trust4.RData")
v.frequency <- data.frame(sample_id = unique(tcr.bcr_trust4$sample_id),
                          IGK = NA, IGL = NA, IGH = NA,
                          TRB = NA, TRA = NA, TRG = NA, TRD = NA)

for (i in 1:nrow(v.frequency)) {
  for (j in colnames(v.frequency[,-1])) {
    v.frequency[i,j] <- sum(
      tcr.bcr_trust4$frequency[tcr.bcr_trust4$sample_id %in% 
                                 v.frequency$sample_id[i] & 
                                 grepl(j, tcr.bcr_trust4$V)])
  }
}

rownames(v.frequency) <- v.frequency$sample_id
v.frequency <- v.frequency[,-1]
save(v.frequency,file="v.frequency.RData")
rm(v.frequency,i,j)

## -- region D--
d.frequency <- data.frame(sample_id = unique(tcr.bcr_trust4$sample_id),
                          IGK = NA, IGL = NA, IGH = NA,
                          TRB = NA, TRA = NA, TRG = NA, TRD = NA)

for (i in 1:nrow(d.frequency)) {
  for (j in colnames(d.frequency[,-1])) {
    d.frequency[i,j] <- sum(
      tcr.bcr_trust4$frequency[tcr.bcr_trust4$sample_id %in% 
                                 d.frequency$sample_id[i] & 
                                 grepl(j, tcr.bcr_trust4$D)])
  }
}

rownames(d.frequency) <- d.frequency$sample_id
d.frequency <- d.frequency[,-1]
save(d.frequency,file="d.frequency.RData")
rm(d.frequency,i,j)

## -- regiao J --
j.frequency <- data.frame(sample_id = unique(tcr.bcr_trust4$sample_id),
                          IGK = NA, IGL = NA, IGH = NA,
                          TRB = NA, TRA = NA, TRG = NA, TRD = NA)

for (i in 1:nrow(j.frequency)) {
  for (j in colnames(j.frequency[,-1])) {
    j.frequency[i,j] <- sum(
      tcr.bcr_trust4$frequency[tcr.bcr_trust4$sample_id %in% 
                                 j.frequency$sample_id[i] & 
                                 grepl(j, tcr.bcr_trust4$J)])
  }
}

rownames(j.frequency) <- j.frequency$sample_id
j.frequency <- j.frequency[,-1]
save(j.frequency,file="j.frequency.RData")
rm(list = ls())
################################################################################
# -- heatmap -- ################################################################
library(ComplexHeatmap)
library(stringr)
library(circlize)
library(RColorBrewer)

## -- metadata -- 
load("v.frequency.RData")
load("d.frequency.RData")
load("j.frequency.RData")

metadata <- data.frame(
  sample_id = rownames(v.frequency),
  tissue = as.factor(c("TP","NT","TP")),
  bcr = rowSums(v.frequency[,1:3]+d.frequency[,1:3]+j.frequency[,1:3]),
  tcr = rowSums(v.frequency[,4:7]+d.frequency[,4:7]+j.frequency[,4:7]))

col.seqGreys <- colorRamp2(breaks = seq(0,10,length.out=9),
                           colors = brewer.pal(9,"Greys"))

col.bin.tissue <- c("NT"="white", "TP"="black")

col.ha <- HeatmapAnnotation(tissue=metadata$tissue,
                            freq_bcr=metadata$tcr,
                            freq_tcr=metadata$bcr,
                            col=list(tissue=col.bin.tissue,
                                     freq_bcr=col.seqGreys,
                                     freq_tcr=col.seqGreys))
## -- construcao do data.v
data.v <- v.frequency %>% as.matrix() %>% t()
all(colnames(data.v)%in% metadata$sample_id)

data.v <- log2(data.v + 1)

## -- transformacao z-score
data.v <- as.matrix(t(scale(t(data.v))))
#teste_z <- (data.v[1,]-mean(data.v[1,]))/sd(data.v[1,])
#all(round(data.v[1,],3) == round(teste_z,3))

## -- ordenacao de data - low e high steroid e total de TCR e BCR
data.v <- data.v[,metadata$sample_id]
data.v <- as.matrix(data.v)

## -- heatmap
### -- dados para argumentos
df.v <- as.data.frame(data.v)
df.v$cell_type <- c(rep("B",3),rep("T",4))

ht.v <- Heatmap(data.v,
                top_annotation = col.ha,
                split = df.v$cell_type,
                show_column_names = FALSE,
                cluster_columns = FALSE,
                cluster_rows = F,
                heatmap_legend_param = list(title="TCR and BCR"),
                col=colorRamp2(breaks = seq(-2,2, length.out=9),
                               colors = rev(brewer.pal(9,"BrBG"))),
                row_title = "Region V",
                row_title_gp = gpar(fontsize=6),
                row_title_side = "left",
                row_names_side = "right",
                row_names_gp = gpar(fontsize=8))

rm(data.v,df.v,v.frequency)

### -- region D --

## -- construcao do data.d
data.d <- d.frequency %>% as.matrix() %>% t()
all(colnames(data.d)%in% metadata$sample_id)

data.d <- log2(data.d + 1)

## -- transformacao z-score
data.d <- data.d[rowSums(data.d)!=0,]
data.d <- as.matrix(t(scale(t(data.d))))
#teste_z <- (data[1,]-mean(data[1,]))/sd(data[1,])
#all(data[1,] == teste_z)

## -- ordenacao de data - low e high steroid e total de TCR e BCR
data.d <- data.d[,metadata$sample_id]
data.d <- as.matrix(data.d)

## -- heatmap
### -- dados para argumentos
df.d <- as.data.frame(data.d)
df.d$cell_type <- c("B","T","T")

ht.d <- Heatmap(data.d,
                split = df.d$cell_type,
                show_column_names = FALSE,
                cluster_columns = FALSE,
                cluster_rows = F,
                show_heatmap_legend = FALSE,
                col=colorRamp2(breaks = seq(-2,2, length.out=9),
                               colors = rev(brewer.pal(9,"BrBG"))),
                row_title = "Region D",
                row_title_gp = gpar(fontsize=6),
                row_title_side = "left",
                row_names_side = "right",
                row_names_gp = gpar(fontsize=8))

ht_list = ht.v %v% ht.d
draw(ht_list, merge_legends=TRUE,annotation_legend_side = "top")

rm(d.frequency,data.d,df.d,ht.d,ht.v)

### -- region J --

## -- construcao do data.j
data.j <- j.frequency %>% as.matrix() %>% t()
all(colnames(data.j)%in% metadata$sample_id)

data.j <- log2(data.j + 1)

## -- transformacao z-score
data.j <- as.matrix(t(scale(t(data.j))))
#teste_z <- (data[1,]-mean(data[1,]))/sd(data[1,])
#all(data[1,] == teste_z)

## -- ordenacao de data - low e high steroid e total de TCR e BCR
data.j <- data.j[,metadata$sample_id]
data.j <- as.matrix(data.j)

## -- heatmap
### -- dados para argumentos
df.j <- as.data.frame(data.j)
df.j$cell_type <- c(rep("B",3),rep("T",4))

ht.j <- Heatmap(data.j,
                split = df.j$cell_type,
                show_column_names = FALSE,
                cluster_columns = FALSE,
                cluster_rows = F,
                show_heatmap_legend = FALSE,
                col=colorRamp2(breaks = seq(-2,2, length.out=9),
                               colors = rev(brewer.pal(9,"BrBG"))),
                row_title = "Region J",
                row_title_gp = gpar(fontsize=6),
                row_title_side = "left",
                row_names_side = "right",
                row_names_gp = gpar(fontsize=8))

ht_list = ht_list %v% ht.j
draw(ht_list, merge_legends=TRUE,annotation_legend_side = "top")

rm(list = ls())

