# lifecourse methylome: contact esther.walton@bristol.ac.uk in case of questions

# step 1: within datasets get basic stats (and cell type composition), also separately for males and females

library(pastecs)

#### 1. set your parameter (you only need to make adjustments in this section) -------------------------------------------------------------------

setwd("./")
cohort = "example_cohort"

tp = "example_time_point" 



# set your methylation file
methylation_file=paste0("./beta.",cohort,".rda")


# set your covariate file and your column header names (i.e., what is your ID column called?, etc)
# covariate file format:
# comma-separated (if not, simply change the "sep" flag in the read.table command in line 41)
# individuals in rows; covariates in columns
# only columns needed: ID, Sex, Age (please remove any cell type estimates)

covariate_file="./samplepheno.csv"
ID = "sentrix"  # your ID column name; has to correspond to the IDs used on your methylation file
age = "agev1"   # your age column name (i.e., "Age" or "age", etc)
sex = "Sex"     # your sex column name (i.e., "Sex" or "sex", "gender", etc) ## IMPORTANT: has to be coded M/F


samplesheet=read.table(covariate_file,sep=",",header = T)
samplesheet<-subset(samplesheet, timepoint == tp)

# merge cell type info
cc=read.table(paste0("cellcounts.",cohort,".txt"), header=T)
samplesheet=merge(samplesheet, cc,by.x=ID,by.y="IID")


#2. load the methylation data -------------------------------------------------------------------
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

meth <- loadRData(methylation_file)

meth <- meth[,colnames(meth) %in% samplesheet[,ID]] #keep the samples that are in covariate file
samplesheet = samplesheet[match(colnames(meth), samplesheet[,ID]),,drop=F]

print(all(colnames(meth)==samplesheet[,ID]))


meth=t(meth)

# 3. get stats and save -------------------------------------------------------------------
meth_stats=stat.desc(meth)
# add IQR
test=apply(meth,2,quantile,na.rm=T)
meth_stats=rbind(meth_stats,test[c(2,4),])

nums <- unlist(lapply(samplesheet, is.numeric))
cov_stats=stat.desc(samplesheet[,nums])
# add IQR
test=apply(samplesheet[,nums],2,quantile,na.rm=T)
cov_stats=rbind(cov_stats,test[c(2,4),])


save(meth_stats,cov_stats,file=paste0("meth_",tp,"_stats_",cohort,".RData"))



