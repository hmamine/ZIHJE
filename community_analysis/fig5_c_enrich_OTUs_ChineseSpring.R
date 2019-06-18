#require packages for the analysis the analysis
pkg=c("plyr", "ggplot2", "phyloseq", "data.table", "reshape2", "GUniFrac", "ape", "phytools", "metagenomeSeq", "VennDiagram")
lapply(pkg, library, character.only = TRUE)

rm(list = ls()[grep("*.*", ls())])

palette.phylum <-c("#EF5656","#47B3DA","#F7A415","#88888a","#2BB065")
list.shape<-c(21,19)

theme_new1 <- theme(
	panel.grid.major = element_blank(), 
	panel.grid.minor = element_blank(),
	axis.title.y=element_blank(),
	axis.title.x=element_blank(),
	legend.position="none"
	)


# log2 fold change
	folch=1
# adjusted p-value
	alpha=0.05
# defining lists
	DT_ad1<-vector( "list" )
	DT_ad2<-vector( "list" )
	p_ad<-vector( "list" )
	list.Tm<-list ( "t8dpi","t4dpi" )
	list.Ts<-list ( "local","adjacent" )

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

	physeq=phyloseq( OTU,TAXA,SD,TREE ) 

#removing samples with less than 600 reads
	physeq.600=prune_samples( sample_sums( physeq ) > 600, physeq )
	 
# subset only the leaf samples
	physeq.leaf=subset_samples( physeq.600, Tissue %in% "leaf" )
	
# subset the data 
	physeq.subset=subset_samples( physeq.leaf, Cultivar %in% "ChSp" )

for( i in list.Ts ) 
{
	print ( paste0( "leaf tissue => ", i ))
	subset.i=c( i )
	physeq.i=subset_samples( physeq.subset, Pos %in% subset.i )

		for( j in list.Tm ) 
		{
			print ( paste0( "time point => ", j ))
			subset.j=c( j )
			physeq.j=subset_samples( physeq.i, Time %in% subset.j )
	
			DT.taxa=data.table( tax, keep.rownames=T, key="rn" )

			#computing differentially abundant OTUs 
			m=as( otu_table( physeq.j ), "matrix" ) + 1L
			t=data.frame( as( tax_table( physeq.j ), "matrix" ) )
			T=AnnotatedDataFrame( t )
			s=as( sample_data( physeq.j ), "data.frame" )
			S=AnnotatedDataFrame( s )
	
			obj=newMRexperiment( m, phenoData=S, featureData=T ) 
			p=cumNormStatFast( obj )	
			objTrim=cumNorm( obj, p=p )
			Treatment = pData( obj )$Treatment
			settings = zigControl( maxit=30, verbose=TRUE )
			dsg1=model.matrix( ~0+Treatment, data=s )

			res1=fitZig( obj=objTrim, mod=dsg1, control=settings, useCSSoffset = TRUE)
			zigFit1=res1$fit
			finalMod1=res1$fit$design

			c.mat1 = makeContrasts ( Treatmentzt	-	Treatmentmock, levels = finalMod1)
			fit1 = contrasts.fit( zigFit1, c.mat1 )
			fit1 = eBayes( fit1 )
			DT_1=fit1$coefficients
			DT_1p=fit1$p.value

			DT.zig<-topTable(fit1, number="inf",  adjust.method="BH" )
			DT.zig$rn<-rownames(DT.zig)
			DT.zig<-data.table(DT.zig, key="rn")

			DT.zig<-DT.zig[,c( 1, 4 ,7 )]	
			DT.taxa=DT.taxa[ DT.zig, ]

			ColPhylum<-c("Actinobacteria","Bacteroidetes","Firmicutes","Proteobacteria")

			DT.taxa$ColPhylum <- ifelse(!DT.taxa$Phylum %in% ColPhylum, "Other", ifelse(DT.taxa$Phylum == "Bacteroidetes", "Bacteroidetes", ifelse(DT.taxa$Phylum == "Firmicutes", "Firmicutes", ifelse (DT.taxa$Phylum == "Actinobacteria", "Actinobacteria", "Proteobacteria"))))
			
			DT.taxa$Shp <- ifelse ( DT.taxa$logFC > 0, "Enriched", "Depleted")
	
			DT.sig<-DT.taxa[DT.taxa$P.Value < alpha ]
#			DT.sig$Abs.logFC<- abs(DT.sig$logFC)

	p.MH=ggplot( DT.sig, aes( x= ColPhylum, y= abs(logFC), color = ColPhylum ) )

	p1=p.MH+geom_boxplot()+geom_jitter( aes(shape=Shp), size=2.5 )+
	theme_bw( )+theme_new1+
	scale_shape_manual(values=list.shape)+
	scale_color_manual(values=palette.phylum)+
	scale_y_continuous( limits=c( 0, 7 ), breaks=c( 1,3,5,7 ) )
	
	print(p1)

	p_ad[[i]][[j]]<- p1
	DT_ad1[[i]][[j]]<- DT.taxa
	DT_ad2[[i]][[j]]<- DT.sig
		}
}
	pl4<-p_ad$local$t4dpi
	pl8<-p_ad$local$t8dpi
	pa4<-p_ad$adjacent$t4dpi
	pa8<-p_ad$adjacent$t8dpi

	gridExtra::grid.arrange(pl4, pl8, pa4, pa8, nrow=4)

