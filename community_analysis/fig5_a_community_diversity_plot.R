# written by M. Amine Hassani - ahassani@bot.uni-kiel.de
#required packages
pkg=c("plyr", "ggplot2", "phyloseq", "data.table", "reshape2", "phytools")

lapply(pkg, library, character.only = TRUE)

rm(list = ls()[grep("*.*", ls())])

alphaInd = c("Shannon", "Observed")

color_palette<-c("#000000","#806600","#803300","666666","#EF5656","#47B3DA","#F7A415","#2BB065")

theme_new <- theme (
	panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(),
	legend.position="none",
	axis.text.x = element_text(angle=90, vjust=1)
	)

# upload and prepare phyloseq objects***
mat=read.table( "otu_table.txt", sep="\t", row.names=1, header=T)
mat=as.matrix(mat)
OTU=otu_table(mat, taxa_are_rows=T) 
tax=read.table("taxonomy.txt", sep="\t", row.names=1, header=1)
tax=as.matrix(tax)
TAXA=tax_table(tax)
sd=read.table("sample_data.txt", sep="\t", row.names=1, header=1)
SD=sample_data(sd)
TREE=read.tree("tree.nwk")

physeq= phyloseq(OTU, TAXA, SD,TREE) 

#removing samples with less than 600 reads
physeq1 = prune_samples(sample_sums(physeq) > 600, physeq)

#subset of leaf and root samples
subset1=c("peat","root")
physeq2=subset_samples(physeq1, !Tissue %in% subset1)

#Rarefication to even depth based on smallest sample size
rf=min(sample_sums(physeq2))
physeq.rrf=rarefy_even_depth(physeq2, rf, replace=TRUE, rngseed = 131)

#Ploting species richness	
	p_nat=plot_richness(physeq.rrf,"Concat","Area" ,measures=alphaInd)
	p_nat$layers <- p_nat$layers[-1]
	
	neworder<-c( "MoT0LChSp","ZtT0LChSp","MoT0LChSpAdj","ZtT0LChSpAdj", 		"MoT4LChSp","ZtT4LChSp","MoT4LChSpAdj","ZtT4LChSpAdj",
	"MoT8LChSp","ZtT8LChSp","MoT8LChSpAdj","ZtT8LChSpAdj",
	"MoT0LObx","ZtT0LObx","MoT0LObxAdj","ZtT0LObxAdj", 
	"MoT4LObx","ZtT4LObx","MoT4LObxAdj","ZtT4LObxAdj",
	"MoT8LObx","ZtT8LObx","MoT8LObxAdj","ZtT8LObxAdj",
	"MoT0RootChSp","ZtT0RootChSp","MoT0RootObx","ZtT0RootObx", 
	"MoT8RootChSp","ZtT8RootChSp","ZtT8RootObx","MoT8RootObx",
	"Peat0dpi","Peat0dpiNOB","Peat8dpi","Peat8dpiNOB" )

	p_nat$data$Concat <- ordered(p_nat$data$Concat, levels= neworder)
	p1=p_nat+geom_boxplot(data=p_nat$data, aes(x=Concat, y=value, color=Area))+
	geom_point(size=4)+theme_bw()+
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
	scale_colour_manual(values=color_palette)

# Plotting Observed OTUs	
subp2 =	ggplot(data=p_nat$data[p_nat$data$variable=="Observed",], aes(x=Concat, y=value, color=Area))+
	geom_boxplot(width=1)+
	geom_point(size=4)+theme_bw()+ theme_new +
	scale_colour_manual(values=color_palette) + 
	facet_wrap(~variable)

# Kruskal-Wallis rank sum test on Observed OTUs
ChSp<-data.table(subp2$data[subp2$data$Cultivar=="ChSp",])

ChSp.Loc<-ChSp[ChSp$Area%in%c("InfArea","mockArea"),]
	
	ChSp.Loc.T0<-ChSp.Loc[ChSp.Loc$Time=="t0dpi",]
	Ob.test.ChSp.Loc.T0<-kruskal.test(data=ChSp.Loc.T0, value ~ Concat)

	ChSp.Loc.T4<-ChSp.Loc[ChSp.Loc$Time=="t4dpi",]
	Ob.test.ChSp.Loc.T4<-kruskal.test(data=ChSp.Loc.T4, value ~ Concat)

	ChSp.Loc.T8<-ChSp.Loc[ChSp.Loc$Time=="t8dpi",]
	Ob.test.ChSp.Loc.T8<-kruskal.test(data=ChSp.Loc.T8, value ~ Concat)

ChSp.Adj<-ChSp[!ChSp$Area%in%c("InfArea","mockArea"),]

	ChSp.Adj.T0<-ChSp.Adj[ChSp.Adj$Time=="t0dpi",]
	Ob.test.ChSp.Adj.T0<-kruskal.test(data=ChSp.Adj.T0, value ~ Concat)

	ChSp.Adj.T4<-ChSp.Adj[ChSp.Adj$Time=="t4dpi",]
	Ob.test.ChSp.Adj.T4<-kruskal.test(data=ChSp.Adj.T4, value ~ Concat)

	ChSp.Adj.T8<-ChSp.Adj[ChSp.Adj$Time=="t8dpi",]
	Ob.test.ChSp.Adj.T8<-kruskal.test(data=ChSp.Adj.T8, value ~ Concat)

Obx<-data.table(subp2$data[!subp2$data$Cultivar=="ChSp",])

Obx.Loc<-Obx[Obx$Area%in%c("InfArea","mockArea"),]

	Obx.Loc.T0<-Obx.Loc[Obx.Loc$Time=="t0dpi",]
	Ob.test.Obx.Loc.T0<-kruskal.test(data=Obx.Loc.T0, value ~ Concat)

	Obx.Loc.T4<-Obx.Loc[Obx.Loc$Time=="t4dpi",]
	Ob.test.Obx.Loc.T4<-kruskal.test(data=Obx.Loc.T4, value ~ Concat)

	Obx.Loc.T8<-Obx.Loc[Obx.Loc$Time=="t8dpi",]
	Ob.test.Obx.Loc.T8<-kruskal.test(data=Obx.Loc.T8, value ~ Concat)

Obx.Adj<-Obx[!Obx$Area%in%c("InfArea","mockArea"),]
	Obx.Adj.T0<-Obx.Adj[Obx.Adj$Time=="t0dpi",]
	Ob.test.Obx.Adj.T0<-kruskal.test(data=Obx.Adj.T0, value ~ Concat)
	Obx.Adj.T4<-Obx.Adj[Obx.Adj$Time=="t4dpi",]
	Ob.test.Obx.Adj.T4<-kruskal.test(data=Obx.Adj.T4, value ~ Concat)
	Obx.Adj.T8<-Obx.Adj[Obx.Adj$Time=="t8dpi",]
	Ob.test.Obx.Adj.T8<-kruskal.test(data=Obx.Adj.T8, value ~ Concat)

# Plotting Shannon indices
subp1 =	ggplot(data=p_nat$data[p_nat$data$variable=="Shannon",], aes(x=Concat, y=value, color=Area))+
	geom_boxplot(width=1)+
	geom_point(size=4)+theme_bw()+ theme_new +
	scale_colour_manual(values=color_palette) + 
	facet_wrap(~variable)

# Kruskal-Wallis rank sum test on Shannon indices
ChSp<-data.table(subp1$data[subp1$data$Cultivar=="ChSp",])

ChSp.Loc<-ChSp[ChSp$Area%in%c("InfArea","mockArea"),]

	ChSp.Loc.T0<-ChSp.Loc[ChSp.Loc$Time=="t0dpi",]
	Sh.test.ChSp.Loc.T0<-kruskal.test(data=ChSp.Loc.T0, value ~ Concat)

	ChSp.Loc.T4<-ChSp.Loc[ChSp.Loc$Time=="t4dpi",]
	Sh.test.ChSp.Loc.T4<-kruskal.test(data=ChSp.Loc.T4, value ~ Concat)

	ChSp.Loc.T8<-ChSp.Loc[ChSp.Loc$Time=="t8dpi",]
	Sh.test.ChSp.Loc.T8<-kruskal.test(data=ChSp.Loc.T8, value ~ Concat)

ChSp.Adj<-ChSp[!ChSp$Area%in%c("InfArea","mockArea"),]

	ChSp.Adj.T0<-ChSp.Adj[ChSp.Adj$Time=="t0dpi",]
	Sh.test.ChSp.Adj.T0<-kruskal.test(data=ChSp.Adj.T0, value ~ Concat)

	ChSp.Adj.T4<-ChSp.Adj[ChSp.Adj$Time=="t4dpi",]
	Sh.test.ChSp.Adj.T4<-kruskal.test(data=ChSp.Adj.T4, value ~ Concat)

	ChSp.Adj.T8<-ChSp.Adj[ChSp.Adj$Time=="t8dpi",]
	Sh.test.ChSp.Adj.T8<-kruskal.test(data=ChSp.Adj.T8, value ~ Concat)

Obx<-data.table(subp1$data[!subp1$data$Cultivar=="ChSp",])

Obx.Loc<-Obx[Obx$Area%in%c("InfArea","mockArea"),]

	Obx.Loc.T0<-Obx.Loc[Obx.Loc$Time=="t0dpi",]
	Sh.test.Obx.Loc.T0<-kruskal.test(data=Obx.Loc.T0, value ~ Concat)

	Obx.Loc.T4<-Obx.Loc[Obx.Loc$Time=="t4dpi",]
	Sh.test.Obx.Loc.T4<-kruskal.test(data=Obx.Loc.T4, value ~ Concat)

	Obx.Loc.T8<-Obx.Loc[Obx.Loc$Time=="t8dpi",]
	Sh.test.Obx.Loc.T8<-kruskal.test(data=Obx.Loc.T8, value ~ Concat)

Obx.Adj<-Obx[!Obx$Area%in%c("InfArea","mockArea"),]

	Obx.Adj.T0<-Obx.Adj[Obx.Adj$Time=="t0dpi",]
	Sh.test.Obx.Adj.T0<-kruskal.test(data=Obx.Adj.T0, value ~ Concat)

	Obx.Adj.T4<-Obx.Adj[Obx.Adj$Time=="t4dpi",]
	Sh.test.Obx.Adj.T4<-kruskal.test(data=Obx.Adj.T4, value ~ Concat)

	Obx.Adj.T8<-Obx.Adj[Obx.Adj$Time=="t8dpi",]
	Sh.test.Obx.Adj.T8<-kruskal.test(data=Obx.Adj.T8, value ~ Concat)

gridExtra::grid.arrange(subp2, subp1, ncol=2)

print("done")
