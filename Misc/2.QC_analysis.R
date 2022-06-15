#!/usr/bin/Rscript

if (!requireNamespace("ggplot2")){
    install.packages("ggplot2")
}
library('ggplot2')

# bam_read_count ----------------------------------------------------------

# ridges_plot -------------------------------------------------------------
if (!requireNamespace("ggridges")){
    install.packages("ggridges")
}
if (!requireNamespace("viridis")){
    install.packages("viridis")
}
library('ggridges')
library('viridis')
library("ggsci")
library("wesanderson")
library(scico)
library(RColorBrewer)
library(extrafont)

files = list.files("../FastQC", "mpileup", full.names = T)

#combine readcounts of different samples
df = data.frame()
for(i in 1:length(files)){
    sample.name = strsplit(files[i], "_")[[1]][1]  
    sample.name = strsplit(sample.name, "/")[[1]][3] 
    tmp = read.table(files[i], stringsAsFactors = F, header = T, na.strings = "")
    tmp = cbind(rep(sample.name,nrow(tmp)),tmp)
    df = rbind(df,tmp)
}

#data clean and calculate alle_freq
names(df)[1:2] = c('Sample', "Segment")

Alle_freq = df$altcount / df$depth
Alle_freq[is.na(Alle_freq)] <- 0

Nt = unlist(lapply(table(df$Sample), function(x) {
    1:x
}))
sample_seq <- seq_along(table(df$Sample))
Y = rep(sample_seq, table(df$Sample))
df = cbind(df, Alle_freq, Nt, Y)
df$log10.depth <- log10(df$depth)
head(df)
y.height <- quantile(df$depth, seq(0, 1, 0.01))[96]
y.log.height <- quantile(df$log10.depth, seq(0, 1, 0.01))[96]

df$Segment_sim <- sapply(df$Segment,function(x){
    tmp <- strsplit(x, "-")[[1]]
    return(tmp[length(tmp)])
})
df$Sample <- as.factor(df$Sample)

#combined Plot
# ggplot(data=df) + 
#     geom_ridgeline_gradient(aes(x=Nt,y=Y*y.height,
#         height=depth,group=Sample,fill=Segment),alpha=0.8)+
#     scale_fill_viridis(discrete = TRUE, direction = -1,alpha=0.8)+
#     scale_y_continuous(breaks = (sample_seq)*y.height, 
#         label = levels(df$Sample)) + 
#     scale_x_continuous(breaks=seq(0,13562,1000)) +
#     ylab('Sample') #+theme_ridges()
# ggsave('../FastQC/Coverage.jpeg',device = 'jpeg',height = 7,width = 9)

library(tidyverse)
df_no_na <- as_tibble(df) %>% filter(!grepl("NA", Sample))
#df_no_na <- as_tibble(df) %>% filter(grepl("NA", Sample))
# sample_name <- read_delim("../data/sample_name.txt", delim = " ", col_names = F)
# df_no_na$Sample <- sapply(as.character(df_no_na$Sample), function(x){
#     tmp <- strsplit(x, "-", fixed = T)[[1]]
#     sample_tmp <- tmp[1]
#     barcode_tmp <- strsplit(tmp[length(tmp)], "WSN")[[1]][2]
#     id_tmp <- paste0(sample_tmp, "-", barcode_tmp, ":")
#     return(sample_name$X2[sample_name$X1 == id_tmp])
# })
# df_no_na$Sample <- factor(df_no_na$Sample , levels = sample_name$X2)

# df_no_na$Ref <- grepl("PR8", df_no_na$Segment)
# df_no_na$Ref <- ifelse(df_no_na$Ref, "PR8", "WSN")
df_no_na$Segment_sim <- sapply(df_no_na$Segment, function(x){tmp = strsplit(x, "-")[[1]]; tmp[length(tmp)]})

# length_wsn = df_no_na %>% group_by(Sample,Segment) %>% summarise(n=n()) %>% filter(Sample == "1-210-2233-Barcode1") %>% filter(grepl("WSN", Segment)) %>% .$n %>% sum()
df_plot = df_no_na #%>% filter(grepl("^1", Sample) & Ref != "WSN")
df_plot$Sample = as.character(df_plot$Sample)
(name_old = unique(df_plot$Sample))
df_plot = df_plot %>% filter(Sample %in% c("1-210-2233-Barcode1", "1-210-2233-Barcode2", "1-210-2233-Barcode3")) # subset samples

name_new = c("210", "217", "2231", "2232", "2233")
df_plot$Sample = factor(df_plot$Sample, levels = name_old, labels = name_new)

#df_plot$Nt_PR8 = df_plot$Nt - length_wsn
df_plot$Segment_sim = factor(df_plot$Segment_sim, levels = c("PB2", "PB1","PA","HA","NP","NA","M","NS"))

p1 <- ggplot(data=df_plot) + 
    geom_col(aes(x = Nt, y = log10(depth), fill =  Segment_sim, color = Segment_sim)) +
    scale_fill_nejm(name="Segment")+
    scale_color_nejm(name="Segment")+
    # facet_grid(cols = vars(Ref), rows = vars(Sample))+
    facet_wrap(vars(Sample), ncol = 1)+
    xlab("Nucleotide")+
    theme(legend.position="bottom")+
    theme(plot.background = element_rect(fill = "transparent",colour = NA), legend.background=element_blank(), legend.key = element_blank(), text = element_text(family = "Times"))
ggsave('../results/Coverage_log.png', plot = p1, device = 'png',height = 12.50,width = 33.02, scale = 0.7, bg = "transparent")


p2 <- ggplot(data=df_no_na) + 
    geom_col(aes(x = Nt-length_wsn, y = depth, fill =  Segment, color = Segment)) +
    scale_fill_viridis(discrete = TRUE)+
    scale_color_viridis(discrete = TRUE)+
    facet_wrap(vars(Sample), ncol = 5)+
    theme(legend.position="bottom")
ggsave('../results/Coverage_raw.jpeg', plot = p2, device = 'jpeg',height = 13.82,width = 17.19, scale = 0.7)

# ggplot(data=df) + 
#     geom_col(aes(x = Nt, y = depth, fill =  Segment, color = Segment)) +
#     scale_fill_viridis(discrete = TRUE)+
#     scale_color_viridis(discrete = TRUE)+
#     facet_wrap(vars(Sample), ncol = 5)
# ggsave('../FastQC/Coverage_2.jpeg',device = 'jpeg',height = 20,width = 9)

# ggplot(data=df) + 
#     geom_ridgeline_gradient(aes(x=Nt,y=Y*y.log.height,
#         height=log10.depth,group=Sample,fill=Segment),alpha=0.8)+
#     scale_fill_viridis(discrete = TRUE, direction = -1,alpha=0.8)+
#     scale_y_continuous(breaks = (sample_seq)*y.log.height, 
#         label = levels(df$Sample)) + 
#     scale_x_continuous(breaks=seq(0,13562,1000)) +
#     ylab('Sample') #+theme_ridges()
# ggsave('../FastQC/Coverage_log10.jpeg',device = 'jpeg',height = 7,width = 9)

# ggplot(data=df) +
#     geom_ridgeline_gradient(aes(x=Nt,y=Y,
#         height=Alle_freq,group=Sample,fill=Segment))+
#     scale_fill_viridis(discrete = TRUE, direction = -1)+
#     scale_y_continuous(breaks = sample_seq,label = levels(df$Sample)) +
#     ylab('Sample')
# ggsave('../FastQC/Alle_freq.jpeg',device = 'jpeg',height = 7,width = 9)

# times = round(df$Alle_freq*1000)+1
# df.density = df[rep(1:nrow(df),times),]
# rownames(df.density)=c()

# ggplot(data=df.density,aes(x=Nt,y=Sample,fill=Segment)) + 
#     geom_density_ridges(alpha=0.8)+
#     scale_fill_viridis(discrete = TRUE, direction = -1)+
#     scale_x_continuous(breaks=seq(0,13562,1000))
# ggsave('../FastQC/Alle_freq_smooth.jpeg',
#     device = 'jpeg',height = 7,width = 9,dpi = 300)


