#require packages for the analysis the analysis
pkg=c("plyr", "ggplot2", "phyloseq", "data.table", "reshape2", "GUniFrac", "ape", "phytools", "metagenomeSeq")
lapply(pkg, library, character.only = TRUE)

rm(list = ls()[grep("*.*", ls())])

color.1<-c( "#47B3DA","#F7A415" )

shape.1=c( 17,19,15 )

shape.2=c( 21,24,17,19 )

theme_new <- theme (
	panel.border = element_rect(),
	panel.grid.major = element_blank(),
	panel.grid.minor = element_blank()
	)

dist1=(	"bray"	)

list.Cv<-list ( "ChSp" , "Obx" )
list.Tm<-list ( "t0dpi","t4dpi","t8dpi" )

res_ad1<-vector( "list" )
res_ad2<-vector( "list" )

print ("ChSp = Chinese Spring, Obx = Obelisk ") 

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

#removing samples with less than 600 reads
	physeq1=prune_samples( sample_sums( physeq ) > 600, physeq )

#normalization of count reads
	otumat=as( otu_table( physeq1 ), "matrix" )
	mp=newMRexperiment( ( otumat ) )
	physeq_norm=phyloseq( otu_table( MRcounts( cumNorm( mp, p=cumNormStat( mp ) ), norm=TRUE, log=TRUE ),taxa_are_rows=TRUE ), TAXA, SD,TREE )

# subset only the leaf samples
	physeq_norm.leaf=subset_samples( physeq_norm, Tissue %in% "leaf" )

# subset to two separates phyloseq objects, local and adjacent tissues
	subset2=c("InfArea","mockArea")

	physeq_norm.subset.Loc=subset_samples( physeq_norm.leaf, Area %in% subset2 )

df <- as.data.frame(matrix(,ncol=1,nrow=4))

df$V1<-c("explained","residuals","total","pvalue")

for( j in list.Cv )
{
	print( j )

	subset2=c( j )

	physeq_cultivar=subset_samples( physeq_norm.subset.Loc, Cultivar %in% subset2 )

	for( i in list.Tm ) 
	{
		print( i ) 

		#subset to one time point 
		subset4=c( i )
		physeq_tmp=subset_samples( physeq_cultivar, Time %in% subset4 )

		#computing Bray-Curtis distances
		dist_BC=distance( physeq_tmp, dist1 )
	
		sd_tab=data.table( as( sample_data( physeq_tmp ), "data.frame" ),keep.rownames=T, key="rn" )
		
		res_ad1[[j]][[i]]<-with( sd_tab, adonis ( dist_BC ~ Treatment ) )

		tmp<-with( sd_tab, adonis ( dist_BC ~ Treatment ) )

		print(tmp)
		
		output<-as.data.frame( tmp$aov.tab$R2 )	

		pvalue<-as.data.frame( tmp$aov.tab[6] )[1,1]

		output <- rbind( output,pvalue ) 

		df <- cbind( df, output ) 
	}
}
  
rownames(df)<-df[,1]
df<-df[-3,-1]
colnames(df)<-c("ChSp_0dpi","ChSp_4dpi","ChSp_8dpi","Obx_0dpi","Obx_4dpi","Obx_8dpi")

	print ("permanova test for adjacent leaf tissue") 
	print(df)

