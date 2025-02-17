<img align="right" src="../imgs/WGSVariantWorkflow.png" width="300" />

# WGS Variant Pipeline @ St. Jude

<img align="right" src="../imgs/GeneralLogo.svg" width="75" />

For more information regarding the PacBio WGS Variant Pipeline, see the links below:

- [App Note](https://www.pacb.com/wp-content/uploads/Application-note-Consolidated-analysis-tools-with-the-PacBio-WGS-Variant-Pipeline.pdf)
- [GitHub](https://github.com/PacificBiosciences/HiFi-human-WGS-WDL)
- [Webinar](https://youtu.be/inxVEjQhntI?si=pqRH2-JpsrjLu6U7)
- [PRISM Presentation](https://youtu.be/fTwRYcBn6w4?si=8KXXw5px--344gSu)

# 1. Prerequisites
If you have not yet run through the prerequisites on the [home page](https://github.com/RhettRautsaw/StJude_PacBio-WDL-tutorial), please do that before continuing. 

# 2. Download and Setup WDL/Resources
## Download PacBio WGS Variant Pipeline (WDL)
```
git clone \
  --depth 1 --branch v2.0.0 --recursive \
  https://github.com/PacificBiosciences/HiFi-human-WGS-WDL.git \
  HiFi-human-WGS-WDL_v2.0.0
```

## Download Reference Resources
A reference genome and additional resources are needed to run several analyses in the WGS Variant Pipeline including alignment, variant calling, filtering, and annotation. We use GRCh38 and the entire reference data bundle (if needed) can be downloaded from [Zenodo](https://zenodo.org/records/13315674). Luckily, these resources have already been downloaded by the Hartwell Center onto St. Jude's HPC and can be easily accessed – no download needed! 

Reference Bundle Directory:
```
/research/rgs01/applications/hpcf/authorized_apps/hartwell/Automation/REF/wdl-humanwgs.v2.0.0.resources
```

> If you're curious how to download these resources yourself, you can follow the instructions in the `resources/ResourceDataPrep.md` file.

### Setup Reference/Tertiary Map Files
For the WDL workflow to locate these reference resources, two tab-delimited map files are supplied as input. An example of these files can be found in the [`resources`](https://github.com/RhettRautsaw/StJude_PacBio-WDL-tutorial/tree/main/WGS_Variant_Pipeline/resources) directory of this repository. 

- [ref_map.tsv](https://github.com/RhettRautsaw/StJude_PacBio-WDL-tutorial/blob/main/WGS_Variant_Pipeline/resources/ref_map.tsv)
- [tertiary_map.tsv](https://github.com/RhettRautsaw/StJude_PacBio-WDL-tutorial/blob/main/WGS_Variant_Pipeline/resources/tertiary_map.tsv)

These files have also been prepared for you and are on St. Jude's HPC – you can use these files directly – no additional prep needed!

Map Files:
```
/research/rgs01/applications/hpcf/authorized_apps/hartwell/Automation/REF/wdl-humanwgs.v2.0.0.resources/ref_map.tsv
/research/rgs01/applications/hpcf/authorized_apps/hartwell/Automation/REF/wdl-humanwgs.v2.0.0.resources/tertiary_map.tsv
```

## Download Example Data (optional)
Finally, we need some data to analyze! Briefly, the workflow expects ~30x WGS coverage.

If you don't have your own data, then you can use [PacBio's GIAB trio dataset](https://www.pacb.com/connect/datasets/). These files have also been downloaded and placed on St. Jude's HPC for this tutorial.

Example Data uBAM Files:
```
# Proband
/research/rgs01/applications/hpcf/authorized_apps/hartwell/Automation/REF/HG002.hifi_reads.bam

# Father
/research/rgs01/applications/hpcf/authorized_apps/hartwell/Automation/REF/HG003.hifi_reads.bam

# Mother
/research/rgs01/applications/hpcf/authorized_apps/hartwell/Automation/REF/HG004.hifi_reads.bam
```

# 3. Setup Input Files
The primary input for miniwdl and the WGS Variant Pipeline is a JSON file with information such as the sample ID, sample sex, HiFi bam file locations, [Human Phenotype Ontology](https://hpo.jax.org/) (if you want to run tertiary analysis), and the location of the reference/tertiary map files. 

If you are running the GIAB trio as an example, then you can use the example files found in the [`resources`](https://github.com/RhettRautsaw/StJude_PacBio-WDL-tutorial/tree/main/WGS_Variant_Pipeline/resources). If you are running your own data, then these files can be used as a template and the sample ID, sex, and HiFi bam files updated.

The WGS Variant Pipeline has two analysis modes:
- Singleton Analysis
- Family Analysis (cohorts of 3-5 related samples)

## Singleton Analysis Input
[singleton.hpc.inputs.json](https://github.com/RhettRautsaw/StJude_PacBio-WDL-tutorial/blob/main/WGS_Variant_Pipeline/resources/singleton.hpc.inputs.json)
```
{
  "humanwgs_singleton.sample_id": "HG002",
  "humanwgs_singleton.sex": "MALE",
  "humanwgs_singleton.hifi_reads": [
    "/research/rgs01/applications/hpcf/authorized_apps/hartwell/Automation/REF/HG002.hifi_reads.bam"
  ],
  "humanwgs_singleton.phenotypes": "HP:0000001,HP:0001518,HP:0002303",
  "humanwgs_singleton.ref_map_file": "/research/rgs01/applications/hpcf/authorized_apps/hartwell/Automation/REF/wdl-humanwgs.v2.0.0.resources/ref_map.tsv",
  "humanwgs_singleton.tertiary_map_file": "/research/rgs01/applications/hpcf/authorized_apps/hartwell/Automation/REF/wdl-humanwgs.v2.0.0.resources/tertiary_map.tsv",
  "humanwgs_singleton.backend": "HPC",
  "humanwgs_singleton.gpu": false,
  #"humanwgs_singleton.gpuType": "String? (optional)",
  "humanwgs_singleton.preemptible": true
}
```

## Family Analysis Input
[family.hpc.inputs.json](https://github.com/RhettRautsaw/StJude_PacBio-WDL-tutorial/blob/main/WGS_Variant_Pipeline/resources/family.hpc.inputs.json)
```
{
  "humanwgs_family.family": {
    "family_id": "GIAB_trio",
    "samples": [
      {
        "sample_id": "HG002",
        "hifi_reads": [
          "/research/rgs01/applications/hpcf/authorized_apps/hartwell/Automation/REF/HG002.hifi_reads.bam"
        ],
        "affected": true,
        "sex": "MALE",
        "father_id": "HG003",
        "mother_id": "HG004"
      },
      {
        "sample_id": "HG003",
        "hifi_reads": [
          "/research/rgs01/applications/hpcf/authorized_apps/hartwell/Automation/REF/HG003.hifi_reads.bam"
        ],
        "affected": true,
        "sex": "MALE"
      },
      {
        "sample_id": "HG004",
        "hifi_reads": [
          "/research/rgs01/applications/hpcf/authorized_apps/hartwell/Automation/REF/HG004.hifi_reads.bam"
        ],
        "affected": true,
        "sex": "FEMALE"
      }
    ]
  },
  "humanwgs_family.phenotypes": "HP:0000001,HP:0001518,HP:0002303",
  "humanwgs_family.ref_map_file": "/research/rgs01/applications/hpcf/authorized_apps/hartwell/Automation/REF/wdl-humanwgs.v2.0.0.resources/ref_map.tsv",
  "humanwgs_family.tertiary_map_file": "/research/rgs01/applications/hpcf/authorized_apps/hartwell/Automation/REF/wdl-humanwgs.v2.0.0.resources/tertiary_map.tsv",
  "humanwgs_family.backend": "HPC",
  "humanwgs_family.gpu": false,
  #"humanwgs_family.gpuType": "String? (optional)",
  "humanwgs_family.preemptible": true
}
```

# 4. Run WGS Variant Pipeline
## Singleton Analysis
miniwdl is a workflow manager that submits a series of parallel jobs to your HPC. Generally, this must be run from the login/head/queue node; therefore, you need to either (1) maintain an active connection until the workflow completes or (2) run the workflow in the background so that if SSH connection to your HPC is lost, the workflow will continue. 

I choose to use `tmux` to setup background jobs; however, St. Jude's HPC allows for sub-job submissions (i.e., jobs to be submitted from compute nodes rather than the head node). Therefore, you are able to use `bsub` instead. 

#### Option 1: bsub
```
bsub -R "rusage[mem=5000]" -q standard \
	miniwdl run HiFi-human-WGS-WDL_v2.0.0/workflows/singleton.wdl \
	--input StJude_PacBio-WDL-tutorial/resources/singleton.hpc.inputs.json \
	--dir pbSingle
```

#### Option 2: tmux
```
tmux new -s pbSingle 

miniwdl run HiFi-human-WGS-WDL_v2.0.0/workflows/singleton.wdl \
	--input StJude_PacBio-WDL-tutorial/resources/singleton.hpc.inputs.json \
	--dir pbSingle
```

Detach tmux session by hitting `Ctrl+b` and then `d`. This will allow miniwdl to continue running in the background. If you'd like to view progress and the output of miniwdl, you can reattach your session by typing:

```
tmux attach -t pbSingle
```

Once miniwdl completes, you can close the tmux session by typing exit in the attached session


## Family Analysis
You can also run the family analysis again either using `tmux` or `bsub`. Below is the base miniwdl code to run.
```
miniwdl run HiFi-human-WGS-WDL_v2.0.0/workflows/family.wdl --input family.hpc.inputs.json --dir pbFamily
```

# 5. Understanding the Output
A directory named `pbSingle`/`pbFamily` will be created and inside these directories will be another dated directory with the format (`YYYYMMDD_HHMMSS_humanwgs_singleton`/`_family`) corresponding to when the workflow was started. If a workflow needs to be restarted, you can submit the same command and it will create a second dated directory and cache the successful parts of the previous run to get to completion faster.

Inside `pbSingle/YYYYMMDD_HHMMSS_humanwgs_singleton`, you will find several `call-*` directories which are the working directories for different parts of the workflow. Unless you are attempting to troubleshoot why your workflow is failing, these can be ignored. 

Focus on the `out` directory and `outputs.json` file as the final outputs. In particular, the `out` directory will contain several sub-directories for different tasks in the workflow. You can find a full description of each of these directories here [(Output Directories Docs)](https://github.com/PacificBiosciences/HiFi-human-WGS-WDL/blob/main/docs/singleton.md).

If this is a family analysis, you may yet find subdirectories beneath this for each sample in your family. These will be numbered based on the order they are supplied in the input WDL. For example:
```
pbFamily/
 └── 20240903_131313_humanwgs_family/
      └── out/
	    └── phased_small_variant_vcf/
		    ├── 0
		    │   └── HG002.GIAB_trio.joint.GRCh38.deepvariant.glnexus.phased.vcf.gz
		    ├── 1
		    │   └── HG003.GIAB_trio.joint.GRCh38.deepvariant.glnexus.phased.vcf.gz
		    └── 2 
			└── HG004.GIAB_trio.joint.GRCh38.deepvariant.glnexus.phased.vcf.gz
```

# 6. Next Steps: Tertiary Analysis
Unfortunately, unlike secondary analysis where we can fairly confidently say that DeepVariant is the recommended tool for SNVs and TRGT is the recommended tool for tandem repeats, tertiary analysis is less straightforward. There are many tools available and no consensus on what tool is best. I would primarily recommend working with our tertiary partners, but there are also a number of third-party software available to test.

This is a rapidly evolving space especially as we learn more and more about the human genome and uncover a more comprehensive view of variants beyond SNVs. This is your space to experiment and see what you can discover!

## PacBio Tertiary Analysis Partners
- [GeneYx](https://geneyx.com/)
- [GoldenHelix](https://www.goldenhelix.com/)

## Third Party Software
### SNVs/Indels
- **Annotation**
	- [slivar](https://github.com/brentp/slivar)
	- [VEP](https://useast.ensembl.org/info/docs/tools/vep/index.html)
	- [SnpEff](https://pcingola.github.io/SnpEff/snpeff/introduction/)/[SnpSift](https://pcingola.github.io/SnpEff/snpsift/introduction/)
	- [oakvar](https://rkimoakbioinformatics.github.io/oakvar/)
- **Pathogenicity prediction**
	- [ClinPred](https://sites.google.com/site/clinpred/)
	- [MAGPIE](https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-023-01274-4)
  
### SVs, CNVs, & Tandem Repeats
- **Annotation**
	- [svpack](https://github.com/PacificBiosciences/svpack)
	- [SnpEff](https://pcingola.github.io/SnpEff/snpeff/introduction/)/[SnpSift](https://pcingola.github.io/SnpEff/snpsift/introduction/)
- **Pathogenicity prediction**
	- [svAnna](https://github.com/monarch-initiative/SvAnna)
	- [Postre](https://github.com/vicsanga/Postre)
	- [TRexs](https://github.com/PacificBiosciences/TRexs) – deprecated
	- [RExPRT](https://github.com/ZuchnerLab/RExPRT)

### Resources/Databases
Useful for filtering out common variants or identifying known pathogenic variants
- [dbSNP](https://www.ncbi.nlm.nih.gov/snp/)
- [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/)
- [dbVar](https://www.ncbi.nlm.nih.gov/dbvar/)
- [OMIM](https://www.omim.org/)
- [gnomAD](https://gnomad.broadinstitute.org/)
- [HPRC](https://humanpangenome.org/data/)
- [AnVIL](https://anvilproject.org/data/consortia?filter=%5B%7B%22categoryKey%22%3A%22dataType%22%2C%22value%22%3A%5B%22Whole+Long-read+Genome%22%2C%22Whole+Long-read+Genome+Methylation%22%5D%7D%5D)
- [colorsDB](https://colorsdb.org/): [Zenodo](https://zenodo.org/records/11511513)
	- [HiFiSolves](https://hifisolves.org/collections)
