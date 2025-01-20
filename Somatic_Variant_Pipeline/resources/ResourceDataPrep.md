# Resource Data Prep

The reference data bundle, input JSON, and example input HiFi data have all already been uploaded to the following directory on St. Jude's HPC.

```
/research/rgs01/applications/hpcf/authorized_apps/hartwell/Automation/REF/wdl-somatic.v0.8.1.resources
```

Therefore, **no further action for data preparation is necessary**!

If there are future releases to the WDL workflow with updated resources, I recommend reaching out to the Hartwell Center to see if they would be able to upload those resources as well. 

Finally, if you would like to download the resources yourself, here are the steps taken:

## Download Reference Resources
Visit [Zenodo](https://doi.org/10.5281/zenodo.8327933) to find the latest version of the reference data bundle and copy the URL. 

You will need to update the URL and name of the resources folder below:
```
cd /research/rgs01/applications/hpcf/authorized_apps/hartwell/Automation/REF/

## Download the reference data bundle
wget https://zenodo.org/records/13347368/files/hifisomatic_resources.tar.gz

## Extract the reference data bundle
tar -xzvf hifisomatic_resources.tar.gz

## Rename the directory and remove the tar file to save room (optional)
mv tarball wdl-somatic.v0.8.1.resources
rm hifisomatic_resources.tar.gz

## Set NEW_RESOURCES variable for next step
NEW_RESOURCES="$PWD/wdl-somatic.v0.8.1.resources"
```

## Download the annotation and hmftools databases
WARNING! These files are very large...only download the example data in a location that has sufficient storage. 

```
cd $NEW_RESOURCES

# Download install script from AnnotSV
wget 'https://github.com/lgmgeo/AnnotSV/raw/master/bin/INSTALL_annotations.sh'
# Download the cache with AnnotSV's script
bash INSTALL_annotations.sh
# The script will create a folder named "AnnotSV_annotations"
# Rename it to just AnnotSV and create a tarball
mv AnnotSV_annotations AnnotSV
tar -czvf annotsv_cache.tar.gz AnnotSV
rm -rf AnnotSV/

# Download VEP bundle
wget 'https://ftp.ensembl.org/pub/release-112/variation/indexed_vep_cache/homo_sapiens_refseq_vep_112_GRCh38.tar.gz'

# Download hmftools resource
wget 'https://storage.googleapis.com/hmf-public/HMFtools-Resources/dna_pipeline/v5_33/38/hmf_dna_pipeline_resources.38_v5.33.tar.gz'

cd ..
```

## Download Example Data (optional)
WARNING! These files are very large...only download the example data in a location that has sufficient storage. 

```
mkdir example_data; cd example_data

# Download tumor demo
wget https://downloads.pacbcloud.com/public/revio/2023Q2/COLO829/COLO829/m84039_230312_025934_s1.hifi_reads.bc2026.bam
ln -s m84039_230312_025934_s1.hifi_reads.bc2026.bam COLO829.tumor.bam

# Download matched normal demo
wget https://downloads.pacbcloud.com/public/revio/2023Q2/COLO829/COLO829-BL/m84039_230327_230708_s1.hifi_reads.bc2007.bam
ln -s m84039_230327_230708_s1.hifi_reads.bc2007.bam COLO829.normal.bam
```
