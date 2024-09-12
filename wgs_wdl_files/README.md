# Resource Data Prep

The reference data bundle, input reference map (`ref_map.tsv`), input tertiary map (`tertiary_map.tsv`), and example input HiFi data have all already been uploaded to the following directory on St. Jude's HPC.

```
/research/rgs01/applications/hpcf/authorized_apps/hartwell/Automation/REF/
```

Therefore, **no further action for data preparation is necessary**!

If there are future releases to the WDL workflow with updated resources, I recommend reaching out to the Hartwell Center to see if they would be able to upload those resources as well. 

Finally, if you would like to download the resources yourself, here are the steps taken:

## Download Reference Resources
Visit [Zenodo](https://zenodo.org/records/13315674) to find the latest version of the reference data bundle and copy the URL. 

You will need to update the URL and name of the resources folder below:
```
## Download the reference data bundle
wget https://zenodo.org/records/13315674/files/hifi-wdl-resources-v2.0.0-rc2.tar

## Extract the reference data bundle
tar -xvf hifi-wdl-resources-v2.0.0-rc2.tar

## Rename the directory and remove the tar file to save room (optional)
mv hifi-wdl-resources-v2.0.0-rc2 wdl-humanwgs.v2.0.0-rc2.resources
rm hifi-wdl-resources-v2.0.0-rc2.tar

## Set NEW_RESOURCES variable for next step
NEW_RESOURCES="$PWD/wdl-humanwgs.v2.0.0-rc2.resources"
```

## Setup Reference/Tertiary Map Files

You can edit the ref_map files yourself or use the files in this folder as a template to replace the old resource path with the new resource path.
```
OLD_RESOURCES="/research/rgs01/applications/hpcf/authorized_apps/hartwell/Automation/REF/wdl-humanwgs.v2.0.0-rc2.resources"

perl -pe "s|$OLD_RESOURCES|$NEW_RESOURCES|g" wgs_wdl_files/ref_map.tsv > $NEW_RESOURCES/ref_map.tsv
perl -pe "s|$OLD_RESOURCES|$NEW_RESOURCES|g" wgs_wdl_files/tertiary_map.tsv > $NEW_RESOURCES/tertiary_map.tsv
```

## Download Example Data (optional)

WARNING! These files are very large...only download the example data in a location that has sufficient storage. 
```
# HG002 (proband)
wget https://downloads.pacbcloud.com/public/revio/2022Q4/HG002-rep1/m84011_220902_175841_s1.hifi_reads.bam
ln -s m84011_220902_175841_s1.hifi_reads.bam HG002.hifi_reads.bam

# HG003 (father)
wget https://downloads.pacbcloud.com/public/revio/2022Q4/HG003-rep1/m84010_220919_235306_s2.hifi_reads.bam
ln -s m84010_220919_235306_s2.hifi_reads.bam HG003.hifi_reads.bam

# HG004 (mother)
wget https://downloads.pacbcloud.com/public/revio/2022Q4/HG004-rep1/m84010_220919_232145_s1.hifi_reads.bam
ln -s m84010_220919_232145_s1.hifi_reads.bam HG004.hifi_reads.bam
```
