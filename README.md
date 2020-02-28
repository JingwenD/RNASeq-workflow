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
22th Feb 2020. Changed the --sjdbGTFfile Homo_sapiens.GRCh38.99.gtf

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
5th Nov. Changed the --outSAMtype as "BAM SortedByCoordinate"

## Feb 22 2020: MAPQ calling via all BAM files

```{bash}
#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -e /hpc/dla_lti/jdeng/logFolder/MAPQcalling.log
#$ -l h_rt=12:00:00
#$ -l h_vmem=1G
#$ -pe threaded 1
#$ -N MAPQcalling
#$ -M jdeng@umcutrecht.nl
#$ -m aes

for i in 103627-001-{81..146}Aligned.out.bam

do

samtools view $i | awk '{print $5}' > $i.csv

done

echo $i
```
# Feb 22 2020: MAPQs sum up

```{bash}
#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -e /hpc/dla_lti/jdeng/logFolder/MAPQstat.log
#$ -l h_rt=02:00:00
#$ -l h_vmem=100G
#$ -pe threaded 2
#$ -N MAPQstat
#$ -M jdeng@umcutrecht.nl
#$ -m aes

cat *.bam.csv | sort -T ./ | uniq -c > MAPQ.txt

```

## Nov 5 2019: Processing BAM with samtools* 
```{bash}
#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -e /hpc/dla_lti/jdeng/logFolder/samtools.log
#$ -l h_rt=24:00:00
#$ -l h_vmem=30G
#$ -pe threaded 6
#$ -N samtools
#$ -M jdeng@umcutrecht.nl
#$ -m aes

for i in *Aligned.out.bam
do
    samtools view -h -q 20 $i |grep -E "^@|NH:i:1$|NH:i:1[^0-9]">$i.unique.sam  #get the read with MAPQ >= 20 and unique mapping from sam file  
    samtools sort -n $i.unique.sam -o $i.sorted.bam   # sorted the sam file with read name and convey to bam file

echo $i
rm $i.unique.sam

done

```


# 3.Counting with *HTSeq* and R package *featureCounts*
Counting assigns mapped sequencing reads to genomic features
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
htseq-count -s no -r name -i gene_id -f bam $i /hpc/dla_lti/jdeng/GenIndex/Homo_sapiens.GRCh38.94.gtf > $i.log

echo $i 

done

```
## Nov 5 2019 Counting with *HTSeq* 
```{bash}
#!/bin/bash 
#$ -S /bin/bash
#$ -cwd
#$ -e /hpc/dla_lti/jdeng/logFolder/htseq2.log
#$ -l h_rt=24:10:00
#$ -l h_vmem=10G
#$ -pe threaded 1
#$ -N htseq2
#$ -M jdeng@umcutrecht.nl
#$ -m aes

module load python/2.7.10

for i in 103627-001-*Aligned.out.bam.sorted.bam; do

htseq-count -s reverse -r name -m intersection-nonempty -a 10 -t exon -i gene_id -f bam $i /hpc/dla_lti/jdeng/GenIndex/Homo_sapiens.GRCh38.94.gtf > $i.countN

echo $i 

done


```

## 3.2 Counting with R package *featureCounts*
```{R}
path <- "/hpc/dla_lti/jdeng/tmp/" 
bam_files <- dir(path) 

counts <- featureCounts("bam_files",

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
write.csv(counts$counts, file = "counts.csv")
write.csv(counts$annotation, file = "annotation.csv")
write.csv(counts$stat, file = "stat.csv")

```

## Arguments
## files	
a character vector giving names of input files containing read mapping results. The files can be in either SAM format or BAM format. The file format is automatically detected by the function.

## annot.inbuilt	
a character string specifying an in-built annotation used for read summarization. It has four possible values including "mm10", "mm9", "hg38" and "hg19", corresponding to the NCBI RefSeq annotations for genomes ‘mm10’, ‘mm9’, ‘hg38’ and ‘hg19’, respectively. "mm10" by default. The in-built annotation has a SAF format (see below).

## annot.ext	
A character string giving name of a user-provided annotation file or a data frame including user-provided annotation data. If the annotation is in GTF format, it can only be provided as a file. If it is in SAF format, it can be provided as a file or a data frame. See below for more details about SAF format annotation. If an annotation file is provided, it can be uncompressed or gzip compressed. Note that annot.ext will override annot.inbuilt if both provided.

## isGTFAnnotationFile	
logical indicating whether the annotation provided via the annot.ext argument is in GTF format or not. FALSE by default. This option is only applicable when annot.ext is not NULL.

## GTF.featureType	
a character string giving the feature type used to select rows in the GTF annotation which will be used for read summarization. "exon" by default. This argument is only applicable when isGTFAnnotationFile is TRUE. Feature types can be found in the third column of a GTF annotation.

## GTF.attrType	
a character string giving the attribute type in the GTF annotation which will be used to group features (eg. exons) into meta-features (eg. genes). "gene_id" by default. This argument is only applicable when isGTFAnnotationFile is TRUE. Attributes can be found in the ninth column of a GTF annotation.

## GTF.attrType.extra	
a character vector specifying extra GTF attribute types that will also be included in the counting output. These attribute types will not be used to group features. NULL by default.

## chrAliases	
a character string giving the name of a chromosome name alias file. This should be a two-column comma-delimited text file. Chromosome name aliases included in this file are used to match chr names in annotation with those in the reads. First column in the file should include chr names in the annotation and second column should include chr names in the reads. Chr names are case sensitive. No column header should be included in the file.

## useMetaFeatures	
logical indicating whether the read summarization should be performed at the feature level (eg. exons) or meta-feature level (eg. genes). If TRUE, features in the annotation (each row is a feature) will be grouped into meta-features, using their values in the GeneID column in the SAF-format annotation file or using the GTF.attrType attribute in the GTF-format annotation file, and then reads will be assiged to the meta-features instead of the features. See below for more details.

## allowMultiOverlap	
logical indicating if a read is allowed to be assigned to more than one feature (or meta-feature) if it is found to overlap with more than one feature (or meta-feature). FALSE by default.

## minOverlap	
integer giving the minimum number of overlapped bases required for assigning a read to a feature (or a meta-feature). For assignment of read pairs (fragments), number of overlapping bases from each read in the same pair will be summed. If a negative value is provided, then a gap of up to specified size will be allowed between read and the feature that the read is assigned to. 1 by default.

## fracOverlap	
numeric giving minimum fraction of overlapping bases in a read that is required for read assignment. Value should be within range [0,1]. 0 by default. Number of overlapping bases is counted from both reads if paired end. Both this option and minOverlap option need to be satisfied before a read can be assigned.

## fracOverlapFeature	
numeric giving minimum fraction of bases included in a feature that is required for overlapping with a read or a read pair. Value should be within range [0,1]. 0 by default.

## largestOverlap	
If TRUE, a read (or read pair) will be assigned to the feature (or meta-feature) that has the largest number of overlapping bases, if the read (or read pair) overlaps with multiple features (or meta-features).

## nonOverlap	
integer giving the maximum number of bases in a read (or a read pair) that are allowed not to overlap the assigned feature. NULL by default (ie. no limit is set).

## nonOverlapFeature	
integer giving the maximum number of non-overlapping bases allowed in a feature during read assignment. NULL by default (ie. no limit is set).

## readShiftType	
a character string indicating how reads should be shifted. It has four possible values including upstream, downstream, left and right. Read shifting is performed before read extension (readExtension5 and readExtension3) and reduction (read2pos).

## readShiftSize	
integer specifying the number of bases the reads will be shifted by. 0 by default. Negative value is not allowed.

## readExtension5	
integer giving the number of bases extended upstream from 5' end of each read. 0 by default. Read extension is performed after read shifting but before read reduction. Negative value is not allowed.

## readExtension3	
integer giving the number of bases extended downstream from 3' end of each read. 0 by default. Negative value is not allowed.

## read2pos	
Specifying whether each read should be reduced to its 5' most base or 3' most base. It has three possible values: NULL, 5 (denoting 5' most base) and 3 (denoting 3' most base). Default value is NULL, ie. no read reduction will be performed. If a read is reduced to a single base, only that base will be considered for the read assignment. Read reduction is performed after read shifting and extension.

## countMultiMappingReads	
logical indicating if multi-mapping reads/fragments should be counted, TRUE by default. ‘NH’ tag is used to located multi-mapping reads in the input BAM/SAM files.

## fraction	
logical indicating if fractional counts are produced for multi-mapping reads and/or multi-overlapping reads. FALSE by default. See below for more details.

## isLongRead	
logical indicating if input data contain long reads. This option should be set to TRUE if counting Nanopore or PacBio long reads.

## minMQS	
integer giving the minimum mapping quality score a read must satisfy in order to be counted. For paired-end reads, at least one end should satisfy this criteria. 0 by default.

## splitOnly	
logical indicating whether only split alignments (their CIGAR strings contain letter 'N') should be included for summarization. FALSE by default. Example split alignments are exon-spanning reads from RNA-seq data. useMetaFeatures should be set to FALSE and allowMultiOverlap should be set to TRUE, if the purpose of summarization is to assign exon-spanning reads to all their overlapping exons.

## nonSplitOnly	
logical indicating whether only non-split alignments (their CIGAR strings do not contain letter 'N') should be included for summarization. FALSE by default.

## primaryOnly	
logical indicating if only primary alignments should be counted. Primary and secondary alignments are identified using bit 0x100 in the Flag field of SAM/BAM files. If TRUE, all primary alignments in a dataset will be counted no matter they are from multi-mapping reads or not (ie. countMultiMappingReads is ignored).

## ignoreDup	
logical indicating whether reads marked as duplicates should be ignored. FALSE by default. Read duplicates are identified using bit Ox400 in the FLAG field in SAM/BAM files. The whole fragment (read pair) will be ignored if paired end.

## strandSpecific	
an integer vector indicating if strand-specific read counting should be performed. Length of the vector should be either 1 (meaning that the value is applied to all input files), or equal to the total number of input files provided. Each vector element should have one of the following three values: 0 (unstranded), 1 (stranded) and 2 (reversely stranded). Default value of this parameter is 0 (ie. unstranded read counting is performed for all input files).

## juncCounts	
logical indicating if number of reads supporting each exon-exon junction will be reported. Junctions are identified from those exon-spanning reads in input data. FALSE by default.

## genome	
a character string giving the name of a FASTA-format file that includes the reference sequences used in read mapping that produced the provided SAM/BAM files. NULL by default. This argument should only be used when juncCounts is TRUE. Note that providing reference sequences is optional when juncCounts is set to TRUE.

## isPairedEnd	
logical indicating if counting should be performed on read pairs or reads. FALSE by default. If TRUE, read pairs will be counted instead of individual reads.

## requireBothEndsMapped	
logical indicating if both ends from the same fragment are required to be successfully aligned before the fragment can be assigned to a feature or meta-feature. This parameter is only appliable when isPairedEnd is TRUE.

## checkFragLength	
logical indicating if the two ends from the same fragment are required to satisify the fragment length criteria before the fragment can be assigned to a feature or meta-feature. This parameter is only appliable when isPairedEnd is TRUE. The fragment length criteria are specified via minFragLength and maxFragLength.

## minFragLength	
integer giving the minimum fragment length for paired-end reads. 50 by default.

## maxFragLength	
integer giving the maximum fragment length for paired-end reads. 600 by default. minFragLength and maxFragLength are only applicable when isPairedEnd is TRUE. Note that when a fragment spans two or more exons, the observed fragment length might be much bigger than the nominal fragment length.

## countChimericFragments	
logical indicating whether a chimeric fragment, which has its two reads mapped to different chromosomes, should be counted or not. TRUE by default.

## autosort	
logical specifying if the automatic read sorting is enabled. This option is only applicable for paired-end reads. If TRUE, reads will be automatically sorted by their names if reads from the same pair are found not to be located next to each other in the input. No read sorting will be performed if there are no such reads found.

## nthreads	
integer giving the number of threads used for running this function. 1 by default.

## byReadGroup	
logical indicating if read counting will be performed for each individual read group. FALSE by default.

## reportReads	
output detailed read assignment results for each read (or fragment if paired end). The detailed assignment results can be saved in three different formats including CORE, SAM and BAM (note that these values are case sensitive). Default value of this option is NULL, indicating not to report assignment results for reads. CORE format represents a tab-delimited text file that contains four columns including read name, status (assigned or the reason if not assigned), number of targets and target list. A target is a feature or a meta-feature. Items in the target lists is separated by comma. If a read is not assigned, its number of targets will be set as -1. For each input file, a text file is generated and its name is the input file name added with ‘.featureCounts’. When SAM or BAM format is specified, the detailed assignment results will be saved to SAM and BAM format files. Names of generated files are the input file names added with ‘.featureCounts.sam’ or ‘.featureCounts.bam’. Three tags are used to describe read assignment results: XS, XN and XT. Tag XS gives the assignment status. Tag XN gives number of targets. Tag XT gives comma separated target list.

## reportReadsPath	
a character string specifying the directory where files including detailed assignment results are saved. If NULL, the results will be saved to the current working directory.

## sampleSheet	
a character string specifying the single-cell RNA sample sheet file. If NULL, featureCounts runs on the bulk RNAseq mode. NULL by default.

## cellBarcodeList	
a character string specifying the file containing the list of cell barcodes for scRNA sample preparation. This argument is mandatory when the sampleSheet option is specified. NULL by default.

## maxMOp	
integer giving the maximum number of ‘M’ operations (matches or mis-matches) allowed in a CIGAR string. 10 by default. Both ‘X’ and ‘=’ operations are treated as ‘M’ and adjacent ‘M’ operations are merged in the CIGAR string.

## tmpDir	
a character string specifying the directory under which intermediate files are saved (later removed). By default, current working directory is used.

## verbose	
logical indicating if verbose information for debugging will be generated. This may include information such as unmatched chromosomes/contigs between reads and annotation.

