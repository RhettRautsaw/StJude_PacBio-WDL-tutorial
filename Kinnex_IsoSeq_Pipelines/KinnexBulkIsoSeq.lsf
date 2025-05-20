#BSUB -J KinnexBulkIsoSeq
#BSUB -q standard
#BSUB -n 64
#BSUB -R "rusage[mem=256000]"
#BSUB -o %J.log

module load parallel/20240222

export HG38=/research/rgs01/applications/hpcf/authorized_apps/hartwell/Automation/REF/Kinnex-IsoSeq/RefGenomes/Human_hg38_Gencode_v39
export REFS=/research/rgs01/applications/hpcf/authorized_apps/hartwell/Automation/REF/Kinnex-IsoSeq

mkdir -p /scratch_space/$USER/KinnexBulkIsoSeqAnalysis
cd /scratch_space/$USER/KinnexBulkIsoSeqAnalysis

singularity pull docker://quay.io/biocontainers/pbskera:1.4.0--hdfd78af_0
singularity pull docker://quay.io/biocontainers/lima:2.13.0--h9ee0642_0
singularity pull docker://quay.io/biocontainers/isoseq:4.3.0--h9ee0642_0
singularity pull docker://quay.io/pacbio/pbmm2:1.17.0_build1
singularity pull docker://quay.io/biocontainers/pbfusion:0.5.1--hdfd78af_0

parallel -a $REFS/bulkRNA/00_example-data/00_samples.txt -j 4 --colsep '\t' '

# Run SKERA
singularity run -B $PWD -B $REFS -B $HG38 docker://quay.io/biocontainers/pbskera:1.4.0--hdfd78af_0 \
	skera split -j 16 $REFS/bulkRNA/00_example-data/heartBulkRNA.kinnex.{1}.bam \
	$REFS/bulkRNA/01_skera-adapters/mas8_primers.fasta \
	01_heartBulkRNA.kinnex.{1}.segmented.bam

# Run LIMA
echo -e "Barcodes,Bio Sample\n{2},{1}" > 02_heartBulkRNA.kinnex.{1}.barcode.csv
singularity run -B $PWD -B $REFS -B $HG38 docker://quay.io/biocontainers/lima:2.13.0--h9ee0642_0 \
	lima --isoseq --log-level INFO -j 16 \
	--overwrite-biosample-names --biosample-csv 02_heartBulkRNA.kinnex.{1}.barcode.csv \
	01_heartBulkRNA.kinnex.{1}.segmented.bam \
	$REFS/bulkRNA/02_lima-primers/IsoSeq_v2_primers_12.fasta \
	02_heartBulkRNA.kinnex.{1}.fl.bam

# Run ISOSEQ REFINE
singularity run -B $PWD -B $REFS -B $HG38 docker://quay.io/biocontainers/isoseq:4.3.0--h9ee0642_0 \
	isoseq refine --require-polya -j 16 \
	02_heartBulkRNA.kinnex.{1}.fl.{2}.bam \
	$REFS/bulkRNA/02_lima-primers/IsoSeq_v2_primers_12.fasta \
	03_heartBulkRNA.kinnex.{1}.flnc.bam

realpath 03_heartBulkRNA.kinnex.{1}*.bam > 03_heartBulkRNA.kinnex.{1}.flnc.fofn

# Run PBMM2
singularity run -B $PWD -B $REFS -B $HG38 docker://quay.io/pacbio/pbmm2:1.17.0_build1 \
	pbmm2 align -j 16 --preset ISOSEQ --sort \
	$HG38/human_GRCh38_no_alt_analysis_set.fasta \
	03_heartBulkRNA.kinnex.{1}.flnc.fofn \
	04_heartBulkRNA.kinnex.{1}.align.bam

# Run pbfusion
singularity run -B $PWD -B $REFS -B $HG38 docker://quay.io/biocontainers/pbfusion:0.5.1--hdfd78af_0 \
	pbfusion discover -t 16 \
	--gtf $HG38/gencode.v39.annotation.sorted.gtf \
	--output-prefix 08_heartBulkRNA.kinnex.{1}.pbfusion \
	04_heartBulkRNA.kinnex.{1}.align.bam
'

realpath 04_heartBulkRNA.kinnex.*.align.bam > 04_heartBulkRNA.kinnex.align.fofn
perl -pe 's/.*kinnex.//g' 04_heartBulkRNA.kinnex.align.fofn | perl -pe 's/.align.bam//g' > 04_heartBulkRNA.kinnex.align.labels

singularity run -B $PWD -B $REFS -B $HG38 docker://quay.io/pacbio/pb_wdl_base:build3 \
	python $REFS/isoquant_generateYAML.py -b 04_heartBulkRNA.kinnex.align.fofn -l 04_heartBulkRNA.kinnex.align.labels -e 05_heartBulkRNA.kinnex.isoquant -o 05_isoquant.yaml

singularity run -B $PWD -B $REFS -B $HG38 docker://quay.io/biocontainers/isoquant:3.6.3--hdfd78af_0 \
	isoquant.py -t 64 -d pacbio --yaml 05_isoquant.yaml \
	-r $HG38/human_GRCh38_no_alt_analysis_set.fasta \
	-g $HG38/gencode.v39.annotation.sorted.gtf.db --complete_genedb \
	-o 05_heartBulkRNA.kinnex.isoquant \
	--sqanti_output

ln -s 05_heartBulkRNA.kinnex.isoquant/05_heartBulkRNA.kinnex.isoquant/05_heartBulkRNA.kinnex.isoquant.* .

singularity run -B $PWD -B $REFS -B $HG38 docker://quay.io/pacbio/pb_wdl_base:build3 \
	python $REFS/isoquant2pigeon.py \
	--gtf 05_heartBulkRNA.kinnex.isoquant.transcript_models.gtf \
	--tsv 05_heartBulkRNA.kinnex.isoquant.transcript_model_grouped_counts.tsv \
	--output 06_heartBulkRNA.kinnex.pigeon.transcript_model_grouped_counts.csv

singularity run -B $PWD -B $REFS -B $HG38 docker://quay.io/biocontainers/pbpigeon:1.4.0--h9948957_0 \
	pigeon classify -j 64 -o 06_heartBulkRNA.kinnex.pigeon \
	05_heartBulkRNA.kinnex.isoquant.transcript_models.gtf \
	$HG38/gencode.v39.annotation.sorted.gtf \
	$HG38/human_GRCh38_no_alt_analysis_set.fasta \
	--flnc 06_heartBulkRNA.kinnex.pigeon.transcript_model_grouped_counts.csv \
	--cage-peak $HG38/refTSS_v3.3_human_coordinate.hg38.sorted.bed \
	--poly-a $HG38/polyA.list.txt \
	--coverage $HG38/intropolis.v1.hg19_with_liftover_to_hg38.tsv.min_count_10.modified2.sorted.tsv

singularity run -B $PWD -B $REFS -B $HG38 docker://quay.io/biocontainers/pbpigeon:1.4.0--h9948957_0 \
	pigeon report -j 64 06_heartBulkRNA.kinnex.pigeon_classification.txt 06_heartBulkRNA.kinnex.pigeon_classification.report.txt

awk '{print $3, $1}' $REFS/bulkRNA/00_example-data/00_samples.txt > 07_samples.matrix.txt
ln -s 05_heartBulkRNA.kinnex.isoquant.transcript_model_grouped_counts.tsv 07_isoquant.isoforms.matrix
ln -s 05_heartBulkRNA.kinnex.isoquant.gene_grouped_counts.tsv 07_isoquant.genes.matrix

singularity run -B $PWD -B $REFS -B $HG38 docker://quay.io/biocontainers/trinity:2.15.2--pl5321h077b44d_3 \
	run_DE_analysis.pl \
	--matrix 07_isoquant.isoforms.matrix \
	--method DESeq2 \
	--samples_file 07_samples.matrix.txt \
	--output 07_deseq2.isoforms

head -n21 07_deseq2.isoforms/*.DE_results | cut -f1 | tail -n +2 > 07_deseq2.top20.isoforms.txt

for i in `cat 07_deseq2.top20.isoforms.txt`
	do 
	grep $i 06_heartBulkRNA.kinnex.pigeon_classification.txt >> 07_deseq2.top20.isoforms.classify.txt
	done

singularity run -B $PWD -B $REFS -B $HG38 docker://quay.io/biocontainers/trinity:2.15.2--pl5321h077b44d_3 \
	run_DE_analysis.pl \
	--matrix 07_isoquant.genes.matrix \
	--method DESeq2 \
	--samples_file 07_samples.matrix.txt \
	--output 07_deseq2.genes

head -n21 07_deseq2.genes/*.DE_results | cut -f1 | tail -n +2 > 07_deseq2.top20.genes.txt

for i in `cat 07_deseq2.top20.genes.txt`
	do 
	grep $i 05_heartBulkRNA.kinnex.isoquant.transcript_models.gtf | grep -P "\tgene\t" >> 07_deseq2.top20.genes.classify.txt
	done
