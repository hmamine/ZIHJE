#require packages
pkg=c("plyr", "ggplot2", "phyloseq", "data.table", "reshape2", "GUniFrac", "ape", "phytools", "metagenomeSeq")
lapply(pkg, library, character.only = TRUE)

rm(list = ls()[grep("*.*", ls())])

color.1<-c( "#47B3DA","#ef5656" )

shape.1=c( 17,19,15 )
shape.2=c( 21,24,17,19 )

theme_new <- theme (
	panel.border = element_rect(),
	panel.grid.major = element_blank(),
	panel.grid.minor = element_blank()
	)

meth1=(	"PCoA"	)
dist1=(	"bray"	)

# upload and prepare phyloseq objects***
mat=read.table( "otu_table.txt", sep="\t", row.names=1, header=T )
mat=as.matrix( mat )
OTU=otu_table( mat, taxa_are_rows=T ) 
tax=read.table( "taxonomy.txt", sep="\t", row.names=1, header=1 )
tax=as.matrix( tax )
TAXA=tax_table( tax )
sd=read.table( "sample_data.txt", sep="\t", row.names=1, header=1 )
SD=sample_data( sd )
TREE=read.tree( "tree.nwk" )

physeq= phyloseq( OTU,TAXA,SD,TREE ) 

# removing samples with less than 600 reads
	physeq1=prune_samples( sample_sums( physeq ) > 600, physeq )

# normalization of count reads
	otumat=as( otu_table( physeq1 ), "matrix" )
	mp=newMRexperiment( ( otumat ) )
	physeq_norm=phyloseq( otu_table( MRcounts( cumNorm( mp, p=cumNormStat( mp ) ), norm=TRUE, log=TRUE ),taxa_are_rows=TRUE ), TAXA, SD,TREE )

# subset only the leaf samples
	physeq_norm.leaf=subset_samples( physeq_norm, Tissue %in% "leaf" )

#  subsampling data to Chinese Spring samples
	physeq_norm.subset=subset_samples( physeq_norm.leaf, Cultivar %in% "ChSp" )


# subset to two separates phyloseq objects, local and adjacent tissues
	subset2=c("InfArea","mockArea")
	physeq_norm.subset.Loc=subset_samples( physeq_norm.subset, Area %in% subset2 )
	physeq_norm.subset.Adj=subset_samples(  physeq_norm.subset, ! Area %in% subset2 )

#PCoA computed on Bray-Curtis distances

# local leaf tissue
	BC.Loc=distance( physeq_norm.subset.Loc, dist1 )
	ord_BC.Loc=ordinate( physeq_norm.subset.Loc, meth1 , BC.Loc )

	p.BC.Loc<-plot_ordination( physeq_norm.subset.Loc, ord_BC.Loc, color="Treatment", shape="Time")
	p.BC.Loc$layers<-p.BC.Loc$layers[-1]

	Loc=p.BC.Loc+geom_point(size=3)+theme_bw()+theme_new+
	scale_shape_manual(values=shape.1)+ 
	scale_color_manual(values=color.1)+
	ggtitle("Chinese Spring local leaf tissue")

# adjacent leaf tissue
	BC.Adj=distance( physeq_norm.subset.Adj, dist1 )
	ord_BC.Adj=ordinate( physeq_norm.subset.Adj, meth1 , BC.Adj )
	
	p.BC.Adj<-plot_ordination( physeq_norm.subset.Adj, ord_BC.Adj, color="Treatment", shape="Time")
	p.BC.Adj$layers<-p.BC.Adj$layers[-1]

	Adj=p.BC.Adj+geom_point(size=3)+theme_bw()+theme_new+
	scale_shape_manual(values=shape.1)+ 
	scale_color_manual(values=color.1)+
	ggtitle("Chinese Spring adjacent leaf tissue")

	gridExtra::grid.arrange(Loc, Adj, ncol=2)

print("done")











