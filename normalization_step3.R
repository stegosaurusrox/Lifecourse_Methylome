# set wd

setwd("./")


# define cohort
cohort="example_cohort"

# define number of pcs
npc_common=20
npc_epic=20
npc_450k=12

# Load meffil and set how many cores to use for parallelization
library(meffil)
options(mc.cores=6)



### step 3: harmonized common or epic probes (done at each site)

# define subset and cohort
for (cgset in c("common","450k")){

    # Generate normalized probe values
    
    # don't forget to set number of PCs here
    if (cgset=="450k"){
      pcs=npc_450k} else {
        pcs=npc_epic
        }
      
    norm.objects=get(load(paste("norm-objects",cgset,cohort,"rda", sep=".")))
    
    load(paste("qcsummary.clean",cgset,cohort,"Robj", sep="."))
    norm.beta <- meffil.normalize.samples(norm.objects, cpglist.remove=qc.summary$bad.cpgs$name)
    save(norm.beta, file=paste("beta",cgset,cohort,"rda", sep = "."))
    
    # Generate normalization report
    pcs <- meffil.methylation.pcs(norm.beta)
    norm.summary <- meffil.normalization.summary(norm.objects, pcs=pcs)
    meffil.normalization.report(norm.summary, output.file=paste("normalization.report",cgset,cohort,"html", sep="."))

}

# merge 450k and common into one dataframe
beta.450k=get(load(paste0("beta.450k.",cohort,".rda")))
beta.common=get(load(paste0("beta.common.",cohort,".rda")))

# reduce 450k to probes not in common
beta.450k=beta.450k[!(row.names(beta.450k) %in% row.names(beta.common)),]
dim(beta.450k)

# merge
beta.450k <- beta.450k[,colnames(beta.450k) %in% colnames(beta.common)] #keep the samples that are in covariate file
beta.common = beta.common[,match(colnames(beta.450k), colnames(beta.common)),drop=F]

print(all(colnames(beta.common)==colnames(beta.450k)))

beta=rbind(beta.common,beta.450k)
save(beta, file=paste("beta",cohort,"rda", sep = "."))


cellcounts<-meffil.estimate.cell.counts.from.betas(beta,cell.type.reference="blood gse35069",verbose=T)
cellcounts<-data.frame(IID=row.names(cellcounts),cellcounts)
write.table(cellcounts,paste("cellcounts",cohort,"txt", sep="."),sep="\t",row.names=F,col.names=T,quote=F)
