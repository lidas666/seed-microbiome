#Seed microbiomes promote Astragalus mongholicus seed germination through pathogen suppression and cellulose degradation

#Da Li1,2,3, Weimin Chen1*, Wen Luo1,4, Haofei Zhang1, Yang Liu1,5, Duntao Shu1* and Gehong Wei1*

###########################αdiversity
library(ggplot2)
library(reshape2)
library(vegan)
library(ape)
library(picante)


alpha <- function(x, tree = NULL, base = exp(1)) {
  est <- estimateR(x)
  Richness <- est[1, ]
  Chao1 <- est[2, ]
  ACE <- est[4, ]
  Shannon <- diversity(x, index = 'shannon', base = base)
  Simpson <- diversity(x, index = 'simpson')   
  Pielou <- Shannon / log(Richness, base)
  goods_coverage <- 1 - rowSums(x == 1) / rowSums(x)
  
  result <- data.frame(Richness, Shannon, Simpson, Pielou, Chao1, ACE, goods_coverage)
  if (!is.null(tree)) {
    PD_whole_tree <- pd(x, tree, include.root = FALSE)[1]
    names(PD_whole_tree) <- 'PD_whole_tree'
    result <- cbind(result, PD_whole_tree)
  }
  result
}




baa<-ggplot(df,aes(x = com, y = Richness,fill = germination)) +
  geom_boxplot(width = .7,show.legend = T,
               position = position_dodge(0.7)) +
  theme_bw(base_size = 16) +
  theme(axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.title.x=element_blank(),
        axis.title.y=element_text(colour='black', size=20,vjust = 1.5),#face = "bold"),
        axis.text.y=element_text(colour='black',size=18),#face = "bold"),
        axis.text.x=element_text(colour = "black",size = 18,#face = "bold",
                                 angle = 45,hjust = 0.5,vjust = 0.5),
        strip.text = element_text(colour = "black",size = 15,face = "bold"),
        legend.position = 'top') +
  scale_fill_manual(values = c('Germinating'="#e31a1c",'Ungerminated'="#1f78b4"),
                    name = '') +
  stat_compare_means(aes(group=germination),
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                      symbols = c("***", "**", "*", "NS")),label = "p.signif",
                     label.y = 1500,size = 8) +
  ylim(1,1600)
baa


#####################beta diversity
dis1 <- vegdist(otu,method = 'bray')
adonis_result <- adonis(dis1~sample*day*Germination, group, permutations = 999)
summary(adonis_result)
otuput <- data.frame(adonis_result$aov.tab, check.names = FALSE, stringsAsFactors = FALSE)
otuput <- cbind(rownames(otuput), otuput)
names(otuput) <- c('', 'Df', 'Sums of squares', 'Mean squares', 'F.Model', 'Variation (R2)', 'Pr (>F)')
otuput 

#####################beta dispersion
mod <- betadisper(d = dis1, group = group$Germination, type = 'centroid')
mod
disper<-data.frame(mod$distances)

####################################Edge
library(edgeR)
targets <-read.csv('w.csv',row.names = 1)
targets <-targets[1:8]
targets <- targets[which(rowSums(targets) > 0), ] 
targets <-as.matrix(targets)
group <- rep(c('Ger', 'Unger'), each = 4)
dgelist <- DGEList(counts = targets, group = group)
keep <- rowSums(cpm(dgelist) > 1 ) >= 2
dgelist <- dgelist[keep, ,keep.lib.sizes = FALSE]
dgelist_norm <- calcNormFactors(dgelist, method = 'TMM')
plotMDS(dgelist_norm, col = rep(c('red', 'blue'), each = 5), dim = c(1, 2))
design <- model.matrix(~group)   
dge <- estimateDisp(dgelist_norm, design, robust = TRUE) 
plotBCV(dge) 


fit <- glmFit(dge, design, robust = TRUE)    
lrt <- glmLRT(fit)   

topTags(lrt)
write.csv(topTags(lrt, n = nrow(dgelist$counts)), 'seed4glmLRT.csv', quote = FALSE)       
dge_de <- decideTestsDGE(lrt, adjust.method = 'fdr', p.value = 0.05) 
summary(dge_de)
write.csv(dge_de,"seed4-otu.csv")
plotMD(lrt, status = dge_de, values = c(1, -1), col = c("#e31a1c","#1f78b4"))   
abline(h = c(-1, 1), col = 'gray', lty = 2)

library(ggplot2)
gene <- read.csv('seed4glmLRT.csv',row.names=1)
gene$label <- c(rownames(gene)[1:5],rep(NA,(nrow(gene)-5)))
gene[which(gene$FDR < 0.05 & gene$logFC >= 1),'sig'] <- 'Up'
gene[which(gene$FDR >= 0.05 | abs(gene$logFC) < 1),'sig'] <- 'No sig'
gene[which(gene$FDR < 0.05 & gene$logFC <= -1),'sig'] <- 'Down'

write.csv(gene,"xxx.csv")
ggplot(BRCA_Match_DEG, aes(x =-logFC, y=-log10(FDR), colour=sig)) +
  geom_point(alpha=0.85, size=1.2) +
  scale_color_manual(values=c("#e31a1c",'gray',"#1f78b4")) +
  xlim(c(-11, 11)) + 
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8)+
  geom_hline(yintercept = -log10(0.05), lty=4,col="black",lwd=0.8) + 
  labs(x="logFC", y="-log10(PValue)") +
  ggtitle("Different ASVs") + 
  theme_bw()+
  theme(axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.title.x=element_text( size=18,vjust = 1.5),
        axis.title.y=element_text( size=18,vjust = 1.5),
        axis.text.y=element_text(size=15),
        axis.text.x=element_text(size = 15,
                                 hjust = 1,vjust = 0.5))


##################################Fig 3 net
otu<- read.csv('bawr.csv', row.name = 1, check.names = FALSE)
otu<-otu[1:16]
otu<- otu[which(rowSums(otu) >0.01), ] 
occor <-corr.test(t(otu),use="pairwise",method="spearman",adjust="fdr",alpha=0.05)

occor.r = occor$r 
occor.p = occor$p 
occor.r[occor.p>0.05|abs(occor.r)<0.6] = 0
diag(occor.r)<-0
occor.r

##################################Fig4C, heatmap link
library(vegan)
library(ggplot2)
library(tidyr)
library(dplyr)
library(patchwork)
dif<-read.csv("difwr_mean.csv",row.names = 1)
ko2<-read.csv("ko2_mean.csv",row.names = 1)
net<-rbind(dif,ko2)
hd1 = t(scale(t(dif)))
gene <- row.names(hd1)
hf1 <- data.frame(gene,hd1)
dt1<- hf1 %>% pivot_longer(
  cols = !gene,
  values_drop_na = FALSE,
  names_to = "groups",
  values_to = "expressions")
head(dt1)
dt1$gene <- factor(dt1$gene,
                   levels = rev(unique(dt1$gene)),
                   ordered = T)
dt1$groups <- factor(dt1$groups,
                     levels = unique(dt1$groups),
                     ordered = T)
dt1 %>% 
  group_by(groups) %>%
  mutate(median_expression = median(expressions, na.rm = TRUE)) %>%
  arrange(groups, desc(median_expression)) -> dt1_clustered

p1 <- ggplot(dt1_clustered, aes(x = groups, y = gene, fill = expressions)) +
  geom_tile(show.legend = TRUE) +
  scale_fill_gradient(low = "#0f86a9", high = "#FC8452", name = "Expression") +
  scale_x_discrete(limits = unique(dt1_clustered$groups)) +  
  scale_y_discrete(position = "right") +
  xlab('') + ylab('') +
  theme(axis.text.x.top = element_text(angle = 90, face = "italic", size = 8, hjust = 1, vjust = 0.5),
        axis.text.y.right = element_text(angle = 0, face = "plain", size = 8, hjust = 0, vjust = 0.5),
        panel.background = element_blank(),
        legend.title = element_text(size = 6),
        legend.position = "left") +
  guides(fill = guide_colourbar(direction = "vertical",
                                title.hjust = 0,
                                title.position = "top",
                                ticks.colour = "white",
                                frame.colour = "red",
                                barheight = 5,
                                barwidth = 0.8)) +
  coord_fixed(ratio = 1, expand = TRUE)

p1

hd2 = t(scale(t(ko2)))
gene <- row.names(hd2)
hf2 <- data.frame(gene,hd2)
dt2<- hf2 %>% pivot_longer(
  cols = !gene,
  values_drop_na = FALSE,
  names_to = "groups",
  values_to = "expressions")
head(dt2)
dt2$gene <- factor(dt2$gene,
                   levels = rev(unique(dt2$gene)),
                   ordered = T)
dt2$groups <- factor(dt2$groups,
                     levels = unique(dt2$groups),
                     ordered = T)
dt2 %>% 
  group_by(groups) %>%
  mutate(median_expression = median(expressions, na.rm = TRUE)) %>%
  arrange(groups, desc(median_expression)) -> dt2_clustered
p2 <- ggplot(dt2_clustered, aes(x = groups, y = gene, fill = expressions)) +
  geom_tile(show.legend = TRUE) +
  scale_fill_gradient(low = "#023f75", high = "#c62d17", name = "Expression") +
  scale_x_discrete(limits = unique(dt2_clustered$groups)) +  
  scale_y_discrete(position = "right") +
  xlab('') + ylab('') +
  theme(axis.text.x.top = element_text(angle = 90, face = "italic", size = 8, hjust = 1, vjust = 0.5),
        axis.text.y.right = element_text(angle = 0, face = "plain", size = 8, hjust = 0, vjust = 0.5),
        panel.background = element_blank(),
        legend.title = element_text(size = 6),
        legend.position = "left") +
  guides(fill = guide_colourbar(direction = "vertical",
                                title.hjust = 0,
                                title.position = "top",
                                ticks.colour = "white",
                                frame.colour = "red",
                                barheight = 5,
                                barwidth = 0.8)) +
  coord_fixed(ratio = 1, expand = TRUE)
p2
###net

dif1<-read.csv("difwr.csv",row.names = 1)
ko21<-read.csv("KO2.csv",row.names = 1)
net<-rbind(dif1,ko21)
library(igraph)
library(dplyr)
library(Hmisc)
b<-rcorr(t(net),type="spearman")
rr<-b$r
pp<-b$P
rr[abs(rr)<0.7]<-0
pp<-p.adjust(pp,method="BH")
pp[pp>=0.01&pp<1]<-0  
pp[pp<0.01]<-1
z<-rr*pp
diag(z)<-0
g<-graph.adjacency(z,weighted=TRUE,mode="undirected")
g <- simplify(g)
g
edg <- data.frame(igraph::as_data_frame(g, what = "both")$edges)
nod <- data.frame(igraph::as_data_frame(g, what = "both")$vertices)
write.csv(edg,"edg_.csv")
write.csv(nod,"node_.csv")
net <- read.csv("net.csv")
h1 <- hf1$gene
sou <- net$Source
y1 <- NULL
for (i in sou) {
  ys <- which(h1 == i)
  y1 <- c(y1, ys)
}
net_y1 <- length(h1) + 1 - y1
h2 <- hf2$gene
tar <- net$Target
y2 <- NULL
for (i in tar) {
  yt <- which(h2 == i)
  y2 <- c(y2, yt)
}
net_y2 <- length(h2) + 1 - y2

x1 <- rep(2, length(tar))
x2 <- rep(3, length(sou))

df <- data.frame(Source = sou, x1 = x1, net_y1 = net_y1, 
                 Target = tar, x2 = x2, net_y2 = net_y2, 
                 PN = net$PN)

n <- length(hf1$gene)
p_mid <- ggplot(df) +
  geom_segment(aes(x = x1, y = net_y1, xend = x2, yend = net_y2, color = PN),
               size =0.5) +
  geom_point(aes(x = x1, y = net_y1), size = 2,
             fill = "orange", color = "tomato", stroke = 1, shape = 21) +
  geom_point(aes(x = x2, y = net_y2), size = 2,
             fill = "orange", color = "tomato", stroke = 1, shape = 21) +
  scale_color_manual(values = c("P" = "#FFBE7A", "N" = "#8ECFC9")) +  
  scale_y_continuous(limits = c(1, n), expand = expansion(add = c(0.5, 0.7))) +
  scale_x_continuous(expand = expansion(0, 0.1)) +
  theme_void()

p_mid

pp<-p1+p_mid+p2
pp



###################################Fig5C
library(ggtreeExtra)
library(ggplot2)
library(reshape2)
library(ggtree)
library(treeio)
library(ggstar)
library(ggnewscale)

set.seed(123)
data<-read.csv("tab2.csv")
data1<-data[1:14]
data1<-data1[-9:-13]
tr <- read.tree("./tree2.nwk")
sc1<-scale(as.numeric(t(data1[2])))
sc2<-scale(as.numeric(t(data1[7])))
data1<-cbind(data1,sc1,sc2)
data2<-data[-1:-14]
ro<-data[1]
data2<-t(scale(t(data2)))
data2<-cbind(ro,data2)
data2<-melt(data2)
sc3<-scale(as.numeric(t(data2[3])))
data2<-cbind(data2,sc3)

p <- ggtree(tr, layout="circular", size=0.4) + 
  geom_treescale(x=1, y=1, fontsize=4, linesize=0.5)+
  geom_tiplab(align = F,
              hjust = 0,
              offset = 3)

p <- ggtree(tr, layout="circular", size=0.2) + 
  geom_treescale(x=1, y=5, fontsize=1.5, linesize=0.5)+
  geom_tiplab(align = F,
              hjust = 0,
              offset = 5,
              size=3)
p

p1 <- p + geom_fruit(data=data1,geom=geom_star,mapping=aes(y=tr.tip.label, fill=Phylum,size=sc1,starshape=Degra),position="identity",starstroke=0.2) +
  
  scale_size_continuous(range=c(0., 1),
                        
                        guide=guide_legend(keywidth=1, keyheight=1, override.aes=list(starshape=15), order=2)) +
  
  scale_fill_manual(values=c("#00B76D","#0D0DBB","#FF5900","#C5199E","#65472F"),
                    
                    guide="none")+
  
  scale_starshape_manual(values=c(1, 15),
                         
                         guide=guide_legend(keywidth=0.5, keyheight=0.5, order=1))



p1 

fix(data2)
p2 <- p1 + 
  new_scale_fill() +
  geom_fruit(
    data=data2,
    geom=geom_tile,
    mapping=aes(y=tr.tip.label, x=variable,alpha=sc3),
    pwidth=0.8,
    axis.params=list(
      axis="x", 
      text.angle=-90, 
      hjust=0 
    )
  ) +
  scale_fill_manual(
    values=c("#DDE5FC", "#BBCBFA"),
    guide=guide_legend(
      keywidth=0.5,
      keyheight=0.5, 
      order=2
    ))
p2



p3 <- p2+new_scale_fill() +
  geom_fruit(
    data=data1,
    geom=geom_bar,
    mapping=aes(y=tr.tip.label, x=GI, fill=Phylum),  
    pwidth=0.6,
    stat="identity",
    orientation="y", 
    axis.params=list(
      axis="x",
      text.angle=-45, 
      hjust=0  
    ),
    grid.params=list() 
  ) + 
  scale_fill_manual(
    values=c("#00B76D","#0D0DBB","#FF5900","#C5199E","#65472F"),
    guide=guide_legend(keywidth=0.5, keyheight=0.5, order=6)
  ) +
  theme(
    legend.background=element_rect(fill=NA),
    legend.title=element_text(size=7), 
    legend.text=element_text(size=6), 
    legend.spacing.y = unit(0.02, "cm")  
  ) 
p3

##################################Fig5d-e
library(ggplot2,car)

library(ggpubr)
aa<-read.csv("tab3.csv",row.names = 1)
aa<-aa[-2:-3]
aa<-aa[1:5]
aa <- na.omit(aa)
dat_plot <- reshape2::melt(aa, id = 'Anti')
p <- ggplot(dat_plot, aes(value, Anti)) +
  geom_point(size=6,col="#e31a1c") +
  facet_wrap(~variable, ncol = 4, scale = 'free') +
  geom_smooth(method = 'lm',col=c("black"))+
  theme_bw()+
  theme(axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.title.x=element_blank(),
        axis.title.y=element_text(colour='black', size=20,vjust = 1.5),
        axis.text.y=element_text(colour='black',size=18),
        axis.text.x=element_text(colour = "black",size = 18,
                                 angle = 0,hjust = 0.5,vjust = 0.5),
        strip.text = element_text(colour = "black",size = 15,face = "bold"),
        legend.position = 'top')
p
env <- c('GR',"GF","GI","SA")
R2_adj <- c()
p_value <- c()

for (i in env) {
  fit_stat <- summary(lm(aa[['Anti']]~aa[[i]])) 
  R2_adj <- c(R2_adj, fit_stat$adj.r.squared)  
  p_value <- c(p_value, fit_stat$coefficients[2,4]) 
}

env_stat <- data.frame(row.names = env, R2_adj, p_value)
env_stat 
env_stat$Anti <- max(aa$Anti) * 0.8
for (i in env) env_stat[i,'value'] <- max(aa[[i]]) * 0.8  
env_stat$variable <- rownames(env_stat)
env_stat$label <- paste('R2.adj =', round(env_stat$R2_adj, 5), '\np =', round(env_stat$p_value, 5))  

env_stat <- env_stat[c('Anti', 'variable', 'value', 'label')]
dat_plot$label <- NA
dat_plot <- rbind(dat_plot, env_stat)  

p1<-p + geom_text(data = dat_plot, aes(label = label), size = 6)
p1

aa<-read.csv("tab3.csv",row.names = 1)
aa<-aa[-1:-2]
aa<-aa[1:5]
aa <- na.omit(aa)

dat_plot <- reshape2::melt(aa, id = 'Degradation')

p <- ggplot(dat_plot, aes(value, Degradation)) +
  geom_point(size=6,col="#e31a1c") +
  facet_wrap(~variable, ncol = 4, scale = 'free') +
  geom_smooth(method = 'lm',col=c("black"))+
  theme_bw()+
  theme(axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.title.x=element_blank(),
        axis.title.y=element_text(colour='black', size=20,vjust = 1.5),
        axis.text.y=element_text(colour='black',size=18),
        axis.text.x=element_text(colour = "black",size = 18,
                                 angle = 0,hjust = 0.5,vjust = 0.5),
        strip.text = element_text(colour = "black",size = 15,face = "bold"),
        legend.position = 'top')
p

env <- c('GR',"GF","GI","SA")
R2_adj <- c()
p_value <- c()

for (i in env) {
  fit_stat <- summary(lm(aa[['Degradation']]~aa[[i]]))  
  R2_adj <- c(R2_adj, fit_stat$adj.r.squared) 
  p_value <- c(p_value, fit_stat$coefficients[2,4])  
}

env_stat <- data.frame(row.names = env, R2_adj, p_value)
env_stat  
env_stat$Degradation <- max(aa$Degradation) * 0.8
for (i in env) env_stat[i,'value'] <- max(aa[[i]]) * 0.8  
env_stat$variable <- rownames(env_stat)
env_stat$label <- paste('R2.adj =', round(env_stat$R2_adj, 5), '\np =', round(env_stat$p_value, 5)) 

env_stat <- env_stat[c('Degradation', 'variable', 'value', 'label')]
dat_plot$label <- NA
dat_plot <- rbind(dat_plot, env_stat) 

p2<-p + geom_text(data = dat_plot, aes(label = label), size = 6)
p2
ggsave("degradation.pdf",p2,width = 12, height = 3)

library(cowplot)
pp<-plot_grid(p1,p2,ncol=1)
pp

############################Fig6a-b
library(reshape2)
bb <- read.csv("soilt.csv",header = TRUE,sep = ",")
bb1 <- melt(bb)
bb$Group <- factor(bb$Group,levels = c("CK","Paenibacillus","Bacillus_1","Bacillus_2"))
bb1$Group <- factor(bb1$Group,levels = c("CK","Paenibacillus","Bacillus_1","Bacillus_2"))
library(multcomp)
bb.sample <- colnames(bb)[2:9]
test.b <- c()
for (i in bb.sample) {
  fit1 <- aov(as.formula(sprintf("%s ~ Group",i)),
              data = bb)
  tuk1<-glht(fit1,linfct=mcp(Group="Tukey"))
  res1 <- cld(tuk1,alpah=0.05)
  test.b <- cbind(test.b,res1$mcletters$Letters)
}
colnames(test.b) <- colnames(bb)[2:9]
test.b <- melt(test.b)
colnames(test.b) <- c("Group","variable","value")
library(tidyverse)
test.b1 <- bb %>% gather(variable,value,-Group) %>% group_by(variable,Group) %>% 
  summarise(Max = max(value))
test.b11 <- dcast(test.b1,Group~variable)
for (i in 2:ncol(test.b11)) {
  test.b11[,i] <- test.b11[,i] + max(test.b11[,i])*0.1
}
test.b11 <- melt(test.b11)
test.b1 <- merge(test.b1,test.b11,by = c("variable","Group"))
test.b2 <- merge(test.b,test.b1,by = c("variable","Group"))

library(ggplot2)
cbbPalette <- c("#e31a1c","#1f78b4","#FF5900", "#00B76D","#0D0DBB","#C5199E",
                "#65472F","#A4C9CC", "#DD7694", "#04376E")
ggplot(bb1,aes(Group,value)) + 
  
  geom_boxplot(aes(fill = Group)) +
  geom_text(data = test.b2,aes(x = Group,y = value.y,label = value.x),
            size = 5,color = "black",fontface = "bold") +
  geom_jitter(aes(fill=Group),position = position_jitter(0.2),shape=21,size=4,color="black")+
  scale_fill_manual(values = cbbPalette) +
  ylab("Germination indices") +
  facet_wrap(.~variable,ncol = 4,scales = "free_y") +
  theme_bw()+
  theme(axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.title.x=element_blank(),
        axis.title.y=element_text(colour='black', size=18,face = "bold",vjust = 1.5),
        axis.text.y=element_text(colour='black',size=15),
        axis.text.x=element_text(colour = "black",size = 15,
                                 angle =45,hjust = 0.5,vjust = 0.5),
        strip.text = element_text(colour = "black",size = 15,face = "bold"),
        legend.position = "none")