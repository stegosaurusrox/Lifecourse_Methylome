# set wd

setwd("./")

# Load meffil and set how many cores to use for parallelization
library(meffil)
options(mc.cores=6)

# ### step 2a: merge all qc object files together ("common" probes)

cohort1="example_cohort1"
cohort2="example_cohort2"
cohort3="example_cohort3"
cohort4="example_cohort4"
cohort5="example_cohort5"
cohort6="example_cohort6"

cgset = "common"

qc.objects1=get(load(paste("./qc.objects.clean",cgset,cohort1,"Robj", sep = ".")))
qc.objects2=get(load(paste("./qc.objects.clean",cgset,cohort2,"Robj", sep = ".")))
qc.objects3=get(load(paste("./qc.objects.clean",cgset,cohort3,"Robj", sep = ".")))
qc.objects4=get(load(paste("./qc.objects.clean",cgset,cohort4,"Robj", sep = ".")))
qc.objects5=get(load(paste("./qc.objects.clean",cgset,cohort5,"Robj", sep = ".")))
qc.objects6=get(load(paste("./qc.objects.clean",cgset,cohort6,"Robj", sep = ".")))

add.cohort.to.qc.objects <- function(qc.objects, cohort) {
  for (i in 1:length(qc.objects))
    qc.objects[[i]]$samplesheet$cohort <- cohort
  qc.objects
}

qc.objects1 <- add.cohort.to.qc.objects(qc.objects1, cohort1)
qc.objects2 <- add.cohort.to.qc.objects(qc.objects2, cohort2)
qc.objects3 <- add.cohort.to.qc.objects(qc.objects3, cohort3)
qc.objects4 <- add.cohort.to.qc.objects(qc.objects4, cohort4)
qc.objects5 <- add.cohort.to.qc.objects(qc.objects5, cohort5)
qc.objects6 <- add.cohort.to.qc.objects(qc.objects6, cohort6)

qc.objects <- c(qc.objects1, qc.objects2, qc.objects3, qc.objects4, qc.objects5, qc.objects6)

y <- meffil.plot.pc.fit(qc.objects)
ggsave(y$plot,filename=paste0("pc-fit-combined.",cgset,".pdf"),height=6,width=6)

# note down number of components
norm.objects <- meffil.normalize.quantiles(qc.objects, random.effects="cohort", number.pcs=20)

norm.objects1 <- norm.objects[names(norm.objects) %in% names(qc.objects1)]
save(norm.objects1,file=paste("norm-objects",cgset,cohort1,"rda", sep="."))

norm.objects2 <- norm.objects[names(norm.objects) %in% names(qc.objects2)]
save(norm.objects2,file=paste("norm-objects",cgset,cohort2,"rda", sep="."))

norm.objects3 <- norm.objects[names(norm.objects) %in% names(qc.objects3)]
save(norm.objects3,file=paste("norm-objects",cgset,cohort3,"rda", sep="."))

norm.objects4 <- norm.objects[names(norm.objects) %in% names(qc.objects4)]
save(norm.objects4,file=paste("norm-objects",cgset,cohort4,"rda", sep="."))

norm.objects5 <- norm.objects[names(norm.objects) %in% names(qc.objects5)]
save(norm.objects5,file=paste("norm-objects",cgset,cohort5,"rda", sep="."))

norm.objects6 <- norm.objects[names(norm.objects) %in% names(qc.objects6)]
save(norm.objects6,file=paste("norm-objects",cgset,cohort6,"rda", sep="."))



# ### step 2b: merge all qc object files together ("epic" probes)
# 
cohort1="example_cohort1"
cohort2="example_cohort2"
cohort3="example_cohort3"
cohort5="example_cohort5"

cgset = "epic"

qc.objects1=get(load(paste("./qc.objects.clean",cgset,cohort1,"Robj", sep = ".")))
qc.objects2=get(load(paste("./qc.objects.clean",cgset,cohort2,"Robj", sep = ".")))
qc.objects3=get(load(paste("./qc.objects.clean",cgset,cohort3,"Robj", sep = ".")))
qc.objects5=get(load(paste("./qc.objects.clean",cgset,cohort5,"Robj", sep = ".")))

qc.objects1 <- add.cohort.to.qc.objects(qc.objects1, cohort1)
qc.objects2 <- add.cohort.to.qc.objects(qc.objects2, cohort2)
qc.objects3 <- add.cohort.to.qc.objects(qc.objects3, cohort3)
qc.objects5 <- add.cohort.to.qc.objects(qc.objects5, cohort5)

qc.objects <- c(qc.objects1, qc.objects2, qc.objects3, qc.objects5)

y <- meffil.plot.pc.fit(qc.objects)
ggsave(y$plot,filename=paste0("pc-fit-combined.",cgset,".pdf"),height=6,width=6)

# note down number of components
norm.objects <- meffil.normalize.quantiles(qc.objects, random.effects="cohort", number.pcs=20)

norm.objects1 <- norm.objects[names(norm.objects) %in% names(qc.objects1)]
save(norm.objects1,file=paste("norm-objects",cgset,cohort1,"rda", sep="."))

norm.objects2 <- norm.objects[names(norm.objects) %in% names(qc.objects2)]
save(norm.objects2,file=paste("norm-objects",cgset,cohort2,"rda", sep="."))

norm.objects3 <- norm.objects[names(norm.objects) %in% names(qc.objects3)]
save(norm.objects3,file=paste("norm-objects",cgset,cohort3,"rda", sep="."))

norm.objects5 <- norm.objects[names(norm.objects) %in% names(qc.objects5)]
save(norm.objects5,file=paste("norm-objects",cgset,cohort5,"rda", sep="."))

# ### step 2c: merge all qc object files together ("450k" probes)
# 
cohort1="example_cohort1"
cohort4="example_cohort4"
cohort6="example_cohort6"

cgset = "450k"

qc.objects1=get(load(paste("./qc.objects.clean",cgset,cohort1,"Robj", sep = ".")))
qc.objects4=get(load(paste("./qc.objects.clean",cgset,cohort4,"Robj", sep = ".")))
qc.objects6=get(load(paste("./qc.objects.clean",cgset,cohort6,"Robj", sep = ".")))

qc.objects1 <- add.cohort.to.qc.objects(qc.objects1, cohort1)
qc.objects4 <- add.cohort.to.qc.objects(qc.objects4, cohort4)
qc.objects6 <- add.cohort.to.qc.objects(qc.objects6, cohort6)

qc.objects <- c(qc.objects1, qc.objects4, qc.objects6)

y <- meffil.plot.pc.fit(qc.objects)
ggsave(y$plot,filename=paste0("pc-fit-combined.",cgset,".pdf"),height=6,width=6)

# note down number of components
norm.objects <- meffil.normalize.quantiles(qc.objects, random.effects="cohort", number.pcs=12) # or use "Slide"

norm.objects1 <- norm.objects[names(norm.objects) %in% names(qc.objects1)]
save(norm.objects1,file=paste("norm-objects",cgset,cohort1,"rda", sep="."))

norm.objects4 <- norm.objects[names(norm.objects) %in% names(qc.objects4)]
save(norm.objects4,file=paste("norm-objects",cgset,cohort4,"rda", sep="."))

norm.objects6 <- norm.objects[names(norm.objects) %in% names(qc.objects6)]
save(norm.objects6,file=paste("norm-objects",cgset,cohort6,"rda", sep="."))


