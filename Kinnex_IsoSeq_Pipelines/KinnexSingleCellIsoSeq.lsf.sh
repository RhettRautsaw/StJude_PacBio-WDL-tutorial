#BSUB -J KinnexSingleCellIsoSeq
#BSUB -q standard
#BSUB -n 64
#BSUB -R "rusage[mem=256000]"
#BSUB -o %J.log

module load parallel/20240222

export HG38=/research/rgs01/applications/hpcf/authorized_apps/hartwell/Automation/REF/Kinnex-IsoSeq/RefGenomes/Human_hg38_Gencode_v39
export REFS=/research/rgs01/applications/hpcf/authorized_apps/hartwell/Automation/REF/Kinnex-IsoSeq

mkdir -p /scratch_space/$USER/KinnexSingleCellIsoSeqAnalysis
cd /scratch_space/$USER/KinnexSingleCellIsoSeqAnalysis

singularity pull docker://quay.io/biocontainers/pbskera:1.4.0--hdfd78af_0
singularity pull docker://quay.io/biocontainers/lima:2.13.0--h9ee0642_0
singularity pull docker://quay.io/biocontainers/isoseq:4.3.0--h9ee0642_0
singularity pull docker://quay.io/pacbio/pbmm2:1.17.0_build1
singularity pull docker://quay.io/biocontainers/pbfusion:0.5.1--hdfd78af_0

parallel -a $REFS/scRNA_10x/00_example-data/00_samples.txt -j 4 --colsep '\t' '

# Run SKERA
singularity run -B $PWD -B $REFS -B $HG38 docker://quay.io/biocontainers/pbskera:1.4.0--hdfd78af_0 \
	skera split -j 32 $REFS/scRNA_10x/00_example-data/{1}.scKinnex.bam \
	$REFS/scRNA_10x/01_skera-adapters/mas16_primers.fasta \
	01_{1}.scKinnex.segmented.bam

# Run LIMA
echo -e "Barcodes,Bio Sample\n{2},{1}" > 02_{1}.scKinnex.barcode.csv
singularity run -B $PWD -B $REFS -B $HG38 docker://quay.io/biocontainers/lima:2.13.0--h9ee0642_0 \
	lima --isoseq --log-level INFO -j 32 \
	--overwrite-biosample-names --biosample-csv 02_{1}.scKinnex.barcode.csv \
	01_{1}.scKinnex.segmented.bam \
	$REFS/scRNA_10x/02_lima-primers/10x_5kit_primers.fasta \
	02_{1}.scKinnex.fl.bam

# Run ISOSEQ TAG
singularity run -B $PWD -B $REFS -B $HG38 docker://quay.io/biocontainers/isoseq:4.3.0--h9ee0642_0 \
	isoseq tag -j 32 --design 16B-10U-13X-T \
	02_{1}.scKinnex.fl.5p--3p.bam \
	03_{1}.scKinnex.flt.bam

# Run ISOSEQ REFINE
singularity run -B $PWD -B $REFS -B $HG38 docker://quay.io/biocontainers/isoseq:4.3.0--h9ee0642_0 \
	isoseq refine --require-polya -j 32 \
	03_{1}.scKinnex.flt.bam \
	$REFS/scRNA_10x/02_lima-primers/10x_5kit_primers.fasta \
	04_{1}.scKinnex.fltnc.bam

# Run ISOSEQ CORRECT
singularity run -B $PWD -B $REFS -B $HG38 docker://quay.io/biocontainers/isoseq:4.3.0--h9ee0642_0 \
	isoseq correct -j 32 --method knee \
	--barcodes $REFS/scRNA_10x/05_isoseq-correct/737K_august_2016.txt.gz \
	04_{1}.scKinnex.fltnc.bam \
	05_{1}.scKinnex.fltnc.bc.bam

# Run ISOSEQ BCSTATS
singularity run -B $PWD -B $REFS -B $HG38 docker://quay.io/biocontainers/isoseq:4.3.0--h9ee0642_0 \
	isoseq bcstats -j 32 --json 05_{1}.scKinnex.fltnc.bc.stats.json \
	--output 05_{1}.scKinnex.fltnc.bc.stats.tsv \
	05_{1}.scKinnex.fltnc.bc.bam

# Run PLOT_KNEES
singularity run -B $PWD -B $REFS -B $HG38 docker://quay.io/pacbio/pb_wdl_base:build3 \
	python $REFS/plot_knees.py --tsv 05_{1}.scKinnex.fltnc.bc.stats.tsv \
	--output 05_{1}.scKinnex.fltnc.bc \
	--estimate_percentile 90

# Run ISOSEQ DEDUP
realpath 05_{1}.scKinnex.fltnc.bc.bam > 05_{1}.scKinnex.fltnc.bc.fofn
singularity run -B $PWD -B $REFS -B $HG38 docker://quay.io/biocontainers/isoseq:4.3.0--h9ee0642_0 \
	isoseq groupdedup -j 32 \
	05_{1}.scKinnex.fltnc.bc.fofn \
	06_{1}.scKinnex.fltnc.bc.dedup.bam
#rm 06_{1}.scKinnex.fltnc.bc.dedup.fasta

# Run addXBtag.py
singularity run -B $PWD -B $REFS -B $HG38 docker://quay.io/pacbio/pb_wdl_base:build3 \
	python $REFS/addXBtag.py 06_{1}.scKinnex.fltnc.bc.dedup.bam {1} 06_{1}.scKinnex.fltnc.bc.dedup.xc.bam

# Run PBMM2
singularity run -B $PWD -B $REFS -B $HG38 docker://quay.io/pacbio/pbmm2:1.17.0_build1 \
	pbmm2 align -j 32 --preset ISOSEQ --sort \
	$HG38/human_GRCh38_no_alt_analysis_set.fasta \
	06_{1}.scKinnex.fltnc.bc.dedup.xc.bam \
	07_{1}.scKinnex.fltnc.bc.dedup.align.bam

# Run pbfusion
singularity run -B $PWD -B $REFS -B $HG38 docker://quay.io/biocontainers/pbfusion:0.5.1--hdfd78af_0 \
	pbfusion discover -t 32 \
	--gtf $HG38/gencode.v39.annotation.sorted.gtf \
	--output-prefix 42_{1}.scKinnex.pbfusion \
	07_{1}.scKinnex.fltnc.bc.dedup.align.bam
'

realpath 07_*.scKinnex.fltnc.bc.dedup.align.bam > 07_scKinnex.fltnc.bc.dedup.align.fofn
perl -pe 's/.*07_//g' 07_scKinnex.fltnc.bc.dedup.align.fofn | perl -pe 's/.scKinnex.fltnc.*//g' > 07_scKinnex.fltnc.bc.dedup.align.labels

singularity run -B $PWD -B $REFS -B $HG38 docker://quay.io/pacbio/pb_wdl_base:build3 \
	python $REFS/isoquant_generateYAML.py -e 08_scKinnex.isoquant -o 08_scKinnex.isoquant.yaml \
	-b 07_scKinnex.fltnc.bc.dedup.align.fofn \
	-l 07_scKinnex.fltnc.bc.dedup.align.labels

singularity run -B $PWD -B $REFS -B $HG38 docker://quay.io/biocontainers/isoquant:3.6.3--hdfd78af_0 \
	isoquant.py -t 64 -d pacbio --yaml 08_scKinnex.isoquant.yaml \
	-r $HG38/human_GRCh38_no_alt_analysis_set.fasta \
	-g $HG38/gencode.v39.annotation.sorted.gtf.db --complete_genedb \
	-o 08_scKinnex.isoquant --read_group tag:XB \
	--sqanti_output --bam_tags CB,XB --counts_format both

ln -s 08_scKinnex.isoquant/08_scKinnex.isoquant/08_scKinnex.isoquant.* .
mkdir 08_scKinnex.isoquant.transcript.mtx 08_scKinnex.isoquant.gene.mtx

singularity run -B $PWD -B $REFS -B $HG38 docker://quay.io/biocontainers/isoquant:3.6.3--hdfd78af_0 \
	python $REFS/convert_grouped_counts.py --output 08_scKinnex.isoquant.transcript.mtx/08_scKinnex.isoquant \
	--input 08_scKinnex.isoquant.transcript_model_grouped_counts_linear.tsv \
	--output_format mtx

singularity run -B $PWD -B $REFS -B $HG38 docker://quay.io/biocontainers/isoquant:3.6.3--hdfd78af_0 \
	python $REFS/convert_grouped_counts.py --output 08_scKinnex.isoquant.gene.mtx/08_scKinnex.isoquant \
	--input 08_scKinnex.isoquant.gene_grouped_counts_linear.tsv \
	--output_format mtx

singularity run -B $PWD -B $REFS -B $HG38 docker://quay.io/pacbio/pb_wdl_base:build3 \
	python $REFS/isoquant2pigeon.py \
	--gtf 08_scKinnex.isoquant.transcript_models.gtf \
	--tsv 08_scKinnex.isoquant.transcript_model_counts.tsv \
	--output 09_scKinnex.pigeon.transcript_model_counts.csv

singularity run -B $PWD -B $REFS -B $HG38 docker://quay.io/biocontainers/pbpigeon:1.4.0--h9948957_0 \
	pigeon classify -j 64 -o 09_scKinnex.pigeon \
	08_scKinnex.isoquant.transcript_models.gtf \
	$HG38/gencode.v39.annotation.sorted.gtf \
	$HG38/human_GRCh38_no_alt_analysis_set.fasta \
	--flnc 09_scKinnex.pigeon.transcript_model_counts.csv \
	--cage-peak $HG38/refTSS_v3.3_human_coordinate.hg38.sorted.bed \
	--poly-a $HG38/polyA.list.txt \
	--coverage $HG38/intropolis.v1.hg19_with_liftover_to_hg38.tsv.min_count_10.modified2.sorted.tsv

singularity run -B $PWD -B $REFS -B $HG38 docker://quay.io/biocontainers/pbpigeon:1.4.0--h9948957_0 \
	pigeon report -j 64 09_scKinnex.pigeon_classification.txt 09_scKinnex.pigeon_classification.report.txt


