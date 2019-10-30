# RNASeq workflow




# 1.QC with fastqc and fastqcr

*fastqc* for individual reads quality control analysis and R package *fastqcr* for the overall report generated.

## 1.1 *fastqc* commands in bash file for qsubmit
```{bash}
#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -e /hpc/dla_lti/jdeng/logFolder/FastQC.log
#$ -l h_rt=12:00:00
#$ -l h_vmem=10G
#$ -pe threaded 10
#$ -N FastQC
#$ -M jdeng@umcutrecht.nl
#$ -m aes


for i in /hpc/dla_lti/APandit/Disease_Cohorts/Tissues/PsA_PsO_skin_RNASeq/*.fastq.gz
do
   fastqc -t 10 $i -o /hpc/dla_lti/jdeng/Psoriasis/RNAseq/QC/fastqc -noextract 
done

```

## 1.2 *fastqcr* commands in R script
```{R}

library(fastqcr)

# Aggregating Multiple FastQC Reports into a Data Frame 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Demo QC directory containing zipped FASTQC reports

qc.dir <- system.file("fastqc", package = "fastqcr")
qc <- qc_aggregate(qc.dir)
qc

# Inspecting QC Problems
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# See which modules failed in the most samples
qc_fails(qc, "module")

# Or, see which samples failed the most
qc_fails(qc, "sample")

# Building Multi QC Reports
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
qc_report(qc.dir, result.file = "multi-qc-report" )

# Building One-Sample QC Reports (+ Interpretation)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
qc.file <- system.file("fastqc", "S1_fastqc.zip", package = "fastqcr")
qc_report(qc.file, result.file = "one-sample-report",
          interpret = TRUE)
```


# 2.Alignment with *STAR*

##  2.1 Generating genome index
```{bash}
#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -e /hpc/dla_lti/jdeng/logFolder/buildingGenome.log
#$ -l h_rt=40:00:00
#$ -l h_vmem=100G
#$ -pe threaded 8
#$ -N buildingGenome
#$ -M jdeng@umcutrecht.nl
#$ -m aes

STAR --runThreadN 8 --runMode genomeGenerate --genomeDir /hpc/dla_lti/jdeng/GenIndex --genomeFastaFiles Homo_sapiens.GRCh38.dna.primary_assembly.fa --sjdbGTFfile Homo_sapiens.GRCh38.94.gtf --sjdbOverhang 100

```

##  2.2 Aligning
```{bash}
#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -e /hpc/dla_lti/jdeng/logFolder/AllAlign1un.log
#$ -l h_rt=12:00:00
#$ -l h_vmem=100G
#$ -pe threaded 8
#$ -N AllAlign1un
#$ -M jdeng@umcutrecht.nl
#$ -m aes

STAR --runThreadN 8 --runMode alignReads --genomeDir /hpc/dla_lti/jdeng/GenIndex --readFilesCommand zcat --readFilesIn /hpc/dla_lti/APandit/Disease_Cohorts/Tissues/PsA_PsO_skin_RNASeq/HK3WHDSXX_103627-001-001_TTACCGAC-CGAATACG_L004_R1.fastq.gz /hpc/dla_lti/APandit/Disease_Cohorts/Tissues/PsA_PsO_skin_RNASeq/HK3WHDSXX_103627-001-001_TTACCGAC-CGAATACG_L004_R2.fastq.gz --outFileNamePrefix 103627-001-001a --outSAMtype BAM Unsorted
echo "103627-001-001a"

STAR --runThreadN 8 --runMode alignReads --genomeDir /hpc/dla_lti/jdeng/GenIndex --readFilesCommand zcat --readFilesIn /hpc/dla_lti/APandit/Disease_Cohorts/Tissues/PsA_PsO_skin_RNASeq/HK3WHDSXX_103627-001-002_TCGTCTGA-GTCCTTGA_L004_R1.fastq.gz,/hpc/dla_lti/APandit/Disease_Cohorts/Tissues/PsA_PsO_skin_RNASeq/HKJY3DSXX_103627-001-002_TCGTCTGA-GTCCTTGA_L004_R1.fastq.gz /hpc/dla_lti/APandit/Disease_Cohorts/Tissues/PsA_PsO_skin_RNASeq/HK3WHDSXX_103627-001-002_TCGTCTGA-GTCCTTGA_L004_R2.fastq.gz,/hpc/dla_lti/APandit/Disease_Cohorts/Tissues/PsA_PsO_skin_RNASeq/HKJY3DSXX_103627-001-002_TCGTCTGA-GTCCTTGA_L004_R2.fastq.gz --outFileNamePrefix 103627-001-002 --outSAMtype BAM Unsorted
echo "103627-001-002a"

STAR --runThreadN 8 --runMode alignReads --genomeDir /hpc/dla_lti/jdeng/GenIndex --readFilesCommand zcat --readFilesIn /hpc/dla_lti/APandit/Disease_Cohorts/Tissues/PsA_PsO_skin_RNASeq/HK3WHDSXX_103627-001-003_TTCCAGGT-CAGTGCTT_L004_R1.fastq.gz /hpc/dla_lti/APandit/Disease_Cohorts/Tissues/PsA_PsO_skin_RNASeq/HK3WHDSXX_103627-001-003_TTCCAGGT-CAGTGCTT_L004_R2.fastq.gz --outFileNamePrefix 103627-001-003 --outSAMtype BAM Unsorted
echo "103627-001-003a"

STAR --runThreadN 8 --runMode alignReads --genomeDir /hpc/dla_lti/jdeng/GenIndex --readFilesCommand zcat --readFilesIn /hpc/dla_lti/APandit/Disease_Cohorts/Tissues/PsA_PsO_skin_RNASeq/HKJY3DSXX_103627-001-004_TACGGTCT-TCCATTGC_L004_R1.fastq.gz,/hpc/dla_lti/APandit/Disease_Cohorts/Tissues/PsA_PsO_skin_RNASeq/HK3WHDSXX_103627-001-004_TACGGTCT-TCCATTGC_L004_R1.fastq.gz /hpc/dla_lti/APandit/Disease_Cohorts/Tissues/PsA_PsO_skin_RNASeq/HKJY3DSXX_103627-001-004_TACGGTCT-TCCATTGC_L004_R2.fastq.gz,/hpc/dla_lti/APandit/Disease_Cohorts/Tissues/PsA_PsO_skin_RNASeq/HK3WHDSXX_103627-001-004_TACGGTCT-TCCATTGC_L004_R2.fastq.gz --outFileNamePrefix 103627-001-004 --outSAMtype BAM Unsorted
echo "103627-001-004a"


STAR --runThreadN 8 --runMode alignReads --genomeDir /hpc/dla_lti/jdeng/GenIndex --readFilesCommand zcat --readFilesIn /hpc/dla_lti/APandit/Disease_Cohorts/Tissues/PsA_PsO_skin_RNASeq/HK3WHDSXX_103627-001-005_AAGACCGT-GTCGATTG_L004_R1.fastq.gz,/hpc/dla_lti/APandit/Disease_Cohorts/Tissues/PsA_PsO_skin_RNASeq/HKJY3DSXX_103627-001-005_AAGACCGT-GTCGATTG_L004_R1.fastq.gz /hpc/dla_lti/APandit/Disease_Cohorts/Tissues/PsA_PsO_skin_RNASeq/HK3WHDSXX_103627-001-005_AAGACCGT-GTCGATTG_L004_R2.fastq.gz,/hpc/dla_lti/APandit/Disease_Cohorts/Tissues/PsA_PsO_skin_RNASeq/HKJY3DSXX_103627-001-005_AAGACCGT-GTCGATTG_L004_R2.fastq.gz --outFileNamePrefix 103627-001-005 --outSAMtype BAM Unsorted
echo "103627-001-005a"

STAR --runThreadN 8 --runMode alignReads --genomeDir /hpc/dla_lti/jdeng/GenIndex --readFilesCommand zcat --readFilesIn /hpc/dla_lti/APandit/Disease_Cohorts/Tissues/PsA_PsO_skin_RNASeq/HK3WHDSXX_103627-001-006_CAGGTTCA-ATAACGCC_L004_R1.fastq.gz /hpc/dla_lti/APandit/Disease_Cohorts/Tissues/PsA_PsO_skin_RNASeq/HK3WHDSXX_103627-001-006_CAGGTTCA-ATAACGCC_L004_R2.fastq.gz --outFileNamePrefix 103627-001-006 --outSAMtype BAM Unsorted
echo "103627-001-006a"

STAR --runThreadN 8 --runMode alignReads --genomeDir /hpc/dla_lti/jdeng/GenIndex --readFilesCommand zcat --readFilesIn /hpc/dla_lti/APandit/Disease_Cohorts/Tissues/PsA_PsO_skin_RNASeq/HK3WHDSXX_103627-001-007_TAGGAGCT-GCCTTAAC_L004_R1.fastq.gz,/hpc/dla_lti/APandit/Disease_Cohorts/Tissues/PsA_PsO_skin_RNASeq/HKJY3DSXX_103627-001-007_TAGGAGCT-GCCTTAAC_L004_R1.fastq.gz /hpc/dla_lti/APandit/Disease_Cohorts/Tissues/PsA_PsO_skin_RNASeq/HK3WHDSXX_103627-001-007_TAGGAGCT-GCCTTAAC_L004_R2.fastq.gz,/hpc/dla_lti/APandit/Disease_Cohorts/Tissues/PsA_PsO_skin_RNASeq/HKJY3DSXX_103627-001-007_TAGGAGCT-GCCTTAAC_L004_R2.fastq.gz --outFileNamePrefix 103627-001-007 --outSAMtype BAM Unsorted
echo "103627-001-007a"

STAR --runThreadN 8 --runMode alignReads --genomeDir /hpc/dla_lti/jdeng/GenIndex --readFilesCommand zcat --readFilesIn /hpc/dla_lti/APandit/Disease_Cohorts/Tissues/PsA_PsO_skin_RNASeq/HK3WHDSXX_103627-001-008_TACTCCAG-GGTATAGG_L004_R1.fastq.gz /hpc/dla_lti/APandit/Disease_Cohorts/Tissues/PsA_PsO_skin_RNASeq/HK3WHDSXX_103627-001-008_TACTCCAG-GGTATAGG_L004_R2.fastq.gz --outFileNamePrefix 103627-001-008 --outSAMtype BAM Unsorted
echo "103627-001-008a"

STAR --runThreadN 8 --runMode alignReads --genomeDir /hpc/dla_lti/jdeng/GenIndex --readFilesCommand zcat --readFilesIn /hpc/dla_lti/APandit/Disease_Cohorts/Tissues/PsA_PsO_skin_RNASeq/HK3WHDSXX_103627-001-009_AGTGACCT-TCTAGGAG_L004_R1.fastq.gz /hpc/dla_lti/APandit/Disease_Cohorts/Tissues/PsA_PsO_skin_RNASeq/HK3WHDSXX_103627-001-009_AGTGACCT-TCTAGGAG_L004_R2.fastq.gz --outFileNamePrefix 103627-001-009 --outSAMtype BAM Unsorted
echo "103627-001-009a"

STAR --runThreadN 8 --runMode alignReads --genomeDir /hpc/dla_lti/jdeng/GenIndex --readFilesCommand zcat --readFilesIn /hpc/dla_lti/APandit/Disease_Cohorts/Tissues/PsA_PsO_skin_RNASeq/HK3WHDSXX_103627-001-010_AGCCTATC-TGCGTAAC_L004_R1.fastq.gz,/hpc/dla_lti/APandit/Disease_Cohorts/Tissues/PsA_PsO_skin_RNASeq/HK3L2DSXX_103627-001-010_AGCCTATC-TGCGTAAC_L004_R1.fastq.gz /hpc/dla_lti/APandit/Disease_Cohorts/Tissues/PsA_PsO_skin_RNASeq/HK3WHDSXX_103627-001-010_AGCCTATC-TGCGTAAC_L004_R2.fastq.gz,/hpc/dla_lti/APandit/Disease_Cohorts/Tissues/PsA_PsO_skin_RNASeq/HK3L2DSXX_103627-001-010_AGCCTATC-TGCGTAAC_L004_R2.fastq.gz --outFileNamePrefix 103627-001-010 --outSAMtype BAM Unsorted
echo "103627-001-010a"



```


# 3.Counting with *HTSeq* and R package *featureCounts*
## 3.1 Counting with *HTSeq* 
```{bash}
#!/bin/bash 
#$ -S /bin/bash
#$ -cwd
#$ -e /hpc/dla_lti/jdeng/logFolder/htsequn.log
#$ -l h_rt=24:10:00
#$ -l h_vmem=100G
#$ -pe threaded 10
#$ -N htsequn
#$ -M jdeng@umcutrecht.nl
#$ -m aes

module load python/2.7.10

for i in *Aligned.out.bam; do
htseq-count -s reverse -r name -i gene_id -f bam $i /hpc/dla_lti/jdeng/GenIndex/Homo_sapiens.GRCh38.94.gtf > $i.count

echo $i 

done

```


## 3.2 Counting with R package *featureCounts*
```{R}
featureCounts("bam_files",

	# annotation
	annot.inbuilt="hg38",
	annot.ext="Homo_sapiens.GRCh38.94.gtf",
	isGTFAnnotationFile=TRUE,
	GTF.featureType="gene",
	GTF.attrType="gene_id",
	GTF.attrType.extra="gene_name",
	chrAliases=NULL,
	
	# level of summarization
	useMetaFeatures=TRUE,
	# overlap between reads and features
        allowMultiOverlap=FALSE,
        minOverlap=1,
        fracOverlap=0,
        fracOverlapFeature=0,
        largestOverlap=FALSE,
        nonOverlap=NULL,
        nonOverlapFeature=NULL,
 
	# Read shift, extension and reduction
        readShiftType="upstream",
        readShiftSize=0,
        readExtension5=0,
        readExtension3=0,
        read2pos=NULL,
	
	# multi-mapping reads
        countMultiMappingReads=FALSE,
  
	# fractional counting
        fraction=TRUE,
  
	# long reads
        isLongRead=FALSE,
 
        # read filtering
        minMQS=0,
        splitOnly=FALSE,
        nonSplitOnly=FALSE,
        primaryOnly=FALSE,
        ignoreDup=FALSE,
  
        # strandness
        strandSpecific=0,
  
	# exon-exon junctions
        juncCounts=FALSE,
        genome=NULL,
   
	# parameters specific to paired end reads
        isPairedEnd=TRUE,
        requireBothEndsMapped=TRUE,
        checkFragLength=TRUE,
        minFragLength=50,
        maxFragLength=600,
        countChimericFragments=FALSE,
        autosort=TRUE,
                           
	# number of CPU threads
        nthreads=4,
        # read group
        byReadGroup=FALSE,
        # report assignment result for each read
        reportReads=NULL,
        reportReadsPath=NULL,
        # miscellaneous
        maxMOp=10,
        tmpDir=".",
        verbose=FALSE)
```
