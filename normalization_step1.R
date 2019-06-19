# set wd

setwd("./")

# define cohort
cohort="example_cohort" 

# Load meffil and set how many cores to use for parallelization
library(meffil)
options(mc.cores=6)

# Generate samplesheet
samplesheet <- meffil.create.samplesheet(path_to_idat_files)

# Or read in samplesheet
#samplesheet <- meffil.read.samplesheet(base=".",pattern="mysamplesheet.csv")


#### qc.object for common and 450k probes (if you use EPIC data, substitute "450k" for "epic" in the next line)

for (cgset in c("common","450k")){
  # Background and dye bias correction, sexprediction, cell counts estimates
  qc.objects <- meffil.qc(samplesheet, cell.type.reference="blood gse35069", verbose=TRUE,
                          featureset = cgset)
  
  # don't use this QC object (contains QC fails)
  #save(qc.objects,file="qc.objects.Robj")
  
  # Generate QC report
  qc.parameters <- meffil.qc.parameters(
    beadnum.samples.threshold             = 0.1,
    detectionp.samples.threshold          = 0.1,
    detectionp.cpgs.threshold             = 0.1, 
    beadnum.cpgs.threshold                = 0.1,
    sex.outlier.sd                        = 5,
    snp.concordance.threshold             = 0.95,
    sample.genotype.concordance.threshold = 0.8
  )
  
  
  qc.summary <- meffil.qc.summary(qc.objects, parameters = qc.parameters)
  meffil.qc.report(qc.summary, output.file=paste("qc-report",cgset,cohort,"html", sep="."))
  
  # Remove outlier samples if necessary
  outlier <- qc.summary$bad.samples
  table(outlier$issue)
  
  qc.objects <- meffil.remove.samples(qc.objects, qc.summary$bad.samples$sample.name)
  save(qc.objects,file=paste("qc.objects.clean",cgset,cohort,"Robj", sep = "."))
  
  # Rerun QC summary on clean dataset
  qc.summary <- meffil.qc.summary(qc.objects,parameters=qc.parameters)
  save(qc.summary, file=paste("qcsummary.clean",cgset,cohort,"Robj", sep = "."))
  meffil.qc.report(qc.summary, paste(output.file="qc-report.clean",cgset,cohort,"html", sep = "."))
  
}

### share the following files:
### - all .html files
### - qc.objects.clean.>cgset<.>cohort<.Robj



