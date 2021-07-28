# Phylogenetic-analysis
# Epidemiology and population structure of *Haemophilus influenzae* causing invasive disease

This project was used to understand the population structure of non-typeable *Haemophilus influenzae* (NTHi).

We used R studio and ggplot2 and ggtree packages to observe the phylogenetic distribution of NTHi genomes and its association with metadata, including a clade-related classification based on the presence or abscence of 17 accessory genes (De Chiara, *et al*. Genome sequencing of disease and carriage isolates of NTHi identifies discrete population structure. Proc Natl Acad Sci 2014;111:5439–5444; and Pinto M, *et al*. Insights into the population structure and pan-genome of *H. influenzae*. Infect Genet Evol 2019;67:126–135).


---
# Figure 4. Assembly-based core-SNP phylogenetic tree and clade distribution of nontypeable *H. influenzae*
Authors: Anna Carrera-Salinas, Aida González-Díaz, Laura Calatayud, Julieta Mercado-Maza,
  Carmen Puig, Dàmaris Berbel, Jordi Càmara, Fe Tubau, Imma Grau, M Ángeles Domínguez,
  Carmen Ardanuy, Sara Martí

---

```{r setup, include=FALSE}
knitr::opts_chunk$set(
	echo = TRUE,
	message = TRUE,
	warning = FALSE
)
```


```{r}
library(ggplot2)
library(ggtree)
library(extrafontdb)
library(extrafont)
loadfonts(device = "win", quiet = TRUE)
```

## 1.Load data
```{r}
##Define the working directory
setwd("C:/Users/anna_/IDIBELL/Gonzalez Diaz, Aida - RECERCA/PROYECTOS/Haemophilus/4. Invasive H. influenzae (2014-2019)/WGS/18_HINF_INVASIVAS/AnalisisClados/R")

##Load the tree file obtained using Parsnp from the Harvest suite (.nwk file)
tree<-read.tree("18_HINF_Invasivas-216-parsnp.nwk")
```

## 2.Tree and metadata
```{r}
##Visualize the tree
p <- ggtree(tree,size=0.2,right=TRUE,branch.length = "none") +
  geom_text(aes(label=node),check_overlap=TRUE)+
  theme(text=element_text(family="Arial Nova Light"))

#Rotate the tree branches
p2<- flip(p,214,226)%>% flip(413,295)%>% flip(144,145)%>% flip(88,302)%>% flip(296,298)
p2$layers[[3]] <- NULL
p3 <- rotate(p2, 294) %>% rotate(320)%>% rotate(348)%>% rotate(335)

#Add metadata
metadata<-read.csv("metadata.csv",
                   header=TRUE, sep = ";")
heatmapData=read.csv("metadata.csv",sep = ";", row.names=1)
rn <- rownames(heatmapData)
heatmapData <- as.data.frame(sapply(heatmapData, as.character))
rownames(heatmapData) <- rn
```

## 3.Construction of the phylogenetic tree, including metadata and bootstrap values
```{r}
#Create the tree with metadata
tree_metadata<-gheatmap(p3, heatmapData, offset = 2,  
                        width=5,
                        colnames_position="top", 
                        colnames_angle=90, 
                        hjust=0,
                        color="lightgrey",
                        font.size=3.5,
                        family="Arial Nova Light") +
  geom_treescale(x=0.05, y=230, offset=2, fontsize = 4,family="Arial Nova Light")+
  geom_tiplab(align=TRUE,linetype='dotted',linesize=.25,color="#525252",size=1)+
  geom_nodelab(aes(label=label,x=branch),size=1,vjust=-.5)+
  geom_point2(aes(subset=node==226), color="#66c2a5", alpha=0.5, size=5)+
  geom_point2(aes(subset=node==303), color="#9e0142", alpha=0.5, size=5)+
  geom_point2(aes(subset=node==336), color="#e6f598", alpha=0.5, size=5)+
  geom_point2(aes(subset=node==359), color="#f46d43", alpha=0.5, size=5)+
  geom_point2(aes(subset=node==348), color="#fee08b", alpha=0.5, size=5)+
  geom_point2(aes(subset=node==145), color="#f46d43", alpha=0.5, size=5)+
  geom_point2(aes(subset=node==144), color="#fee08b", alpha=0.5, size=5)+
  geom_point2(aes(subset=node==300), color="#5e4fa2", alpha=0.5, size=5)+
  geom_point2(aes(subset=node==317), color="#5e4fa2", alpha=0.5, size=5)+
  geom_point2(aes(subset=node==380), color="#5e4fa2", alpha=0.5, size=5)+
  geom_point2(aes(subset=node==296), color="#5e4fa2", alpha=0.5, size=5)+
  geom_point2(aes(subset=node==321), color="#5e4fa2", alpha=0.5, size=5)+
  geom_point2(aes(subset=node==413), color="#5e4fa2", alpha=0.5, size=5)+
  geom_point2(aes(subset=node==214), color="#5e4fa2", alpha=0.5, size=5)+
  scale_y_continuous(expand=c(0, 20))+
  scale_fill_manual(name="",values = c("#BA6261", "#65B0CB", "#A1D6DD","#9e0142", "#f46d43", "#fee08b","#e6f598","#66c2a5", "#5e4fa2","red","grey", "white"),limits=c("Carrera-Salinas", "De Chiara", "Pinto","I","II","III","IV","V","VI","Unclassified","Presence","Abscence"))+
  theme(legend.position="none")


```



## Figure 4.Assembly-based core-SNP phylogenetic tree and clade distribution of nontypeable *H. influenzae*

```{r pressure, echo=FALSE}
plot(tree_metadata)
```

