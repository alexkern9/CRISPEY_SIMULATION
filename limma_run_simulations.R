# load packages
library("ggplot2")
library("edgeR")
library("limma")

# parameters
input_samples_info_df_filename = './six_tp_simulation_limma_sample.tsv'
input_counts_df_filename = '.six_tp_simulation_limma_counts.tsv'
OUTPUT_FILENAME = './six_tp_simulation_limma_results.tsv'

READS_NUM_CUTOFF = 500
print("input files:")
print(input_samples_info_df_filename)
print(input_counts_df_filename)

samples <- read.table(input_samples_info_df_filename , header=TRUE, na.strings = "", sep='\t')
rownames(samples) = samples$sampID
samples = samples[ , !(names(samples) %in% c("sampID","X"))]
samples$replicate<-as.factor(samples$replicate)

# loading oligo counts
sample_counts <- read.csv(input_counts_df_filename , header=TRUE, na.strings = "", sep='\t')
rownames(sample_counts) = sample_counts$strain_id
sample_counts = sample_counts[ , !(names(sample_counts) %in% c('barseq',
                                                               'strain_id',
                                                               'agilent_oligo_id',
                                                               'count',
                                                               'oligo_id',
                                                               'SCD_T0_123_A_counts',
                                                               'SCD_T0_123_B_counts',
                                                               'SCD_T0_123_C_counts',
                                                               'SCD_T0_456_A_counts',
                                                               'SCD_T0_456_B_counts',
                                                               'SCD_T0_456_C_counts',
                                                               't0_total_counts'))]

##########################################################################################################################
# limma voom using all time points
##########################################################################################################################

dge <- DGEList(counts=sample_counts)

# filtering oligos with low read number
isexpr <- (rowSums(dge$counts) >= READS_NUM_CUTOFF)
dge <- dge[isexpr , keep.lib.size = FALSE]

dge <- calcNormFactors(dge) # TMM

design <- model.matrix(~generations+ replicate, samples)
#colnames(design)

# counts -> logCPM
v <- voom(dge, design, plot=FALSE)


fit<-lmFit(v,design)
fit<-eBayes(fit)

# this shows what number of effect are differnt from zero for each variable 
# make sure the replicates has many less significant effects than the time variable
print(summary(decideTests(fit)))

# For significant genes for effect A:
res_df_limma = topTable(fit,coef=2, number = 1000000, confint=TRUE)


write.table(res_df_limma, 
            file=paste(OUTPUT_FILENAME, sep=""), 
            sep="\t", quote = FALSE,na = "", row.names=TRUE, col.names = TRUE)

print('Done running limma-voom for: ')
print(OUTPUT_FILENAME)


