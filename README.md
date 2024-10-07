# PacBio WGS Variant Pipeline @ St. Jude

**Rhett Rautsaw**, Field Applications Bioinformatic Support (FABS) Scientist II, PacBio\
**Daniel Darnell**, Sr. Bioinformatics Analyst, Hartwell Center for Biotechnology at St. Jude Children's Research Hospital

⚠️ **THIS DOCUMENT IS NOT YET COMPLETE** ⚠️

Please feel free to take a look and get started by creating your St. Jude HPC account and setting up your environment below.
However, we are still developing the workshop...so stay tuned and look forward to our event on November 12th for more details!

# Summary
This document is designed to guide researchers at St. Jude on how to setup, use, and understand the output of PacBio's WGS Variant Pipeline specifically on St. Jude's HPC system. 

For more information regarding the pipeline or St. Jude's HPC system, see the links below:
- PacBio WGS Variant Pipeline
    - [App Note](https://www.pacb.com/wp-content/uploads/Application-note-Consolidated-analysis-tools-with-the-PacBio-WGS-Variant-Pipeline.pdf)
    - [GitHub](https://github.com/PacificBiosciences/HiFi-human-WGS-WDL)
    - [Webinar](https://youtu.be/inxVEjQhntI?si=pqRH2-JpsrjLu6U7)
    - [PRISM Presentation](https://youtu.be/fTwRYcBn6w4?si=8KXXw5px--344gSu)
- St. Jude HPC
    - [Documentation](https://wiki.stjude.org/display/HPCF/)

If you do not have command line experience or want a more push-button solution. I recommend checking out our PacBio secondary analysis partners at [DNAstack](https://omics.ai/workflows/pacbio/), [DNAnexus](https://www.pacb.com/wp-content/uploads/PacBio-DNAnexus.pdf), or [FormBio](https://www.pacb.com/wp-content/uploads/FORM-Bio-flyer.pdf). 

# 1. Access St. Jude's HPC
## Create an Account
If you do not already have an account for St. Jude's HPC, please visit [ServiceNow](https://stjude.service-now.com/sp) to request a login account.

## Login, start an interactive session, and load Anaconda
```
# Login to St. Jude HPC
ssh username@hpc.stjude.org

# Start an Interactive Session
hpcf_interactive

# Load Anaconda Module
module load conda3/202311

# Prepare Conda Environment
#conda init
conda create -n miniwdl_env pip
conda activate miniwdl_env
```

# 2. Download and Install Dependencies
You will only need to run through this section once! The next time you login to St. Jude's HPC the software does not need to be redownloaded, installed, or setup. 

If you have already downloaded and installed the necessary dependencies, you can skip to [3. Setup Input Files](#3-setup-input-files).

## Install miniwdl and miniwdl-lsf extension
```
pip install miniwdl

pip install git+https://github.com/adthrasher/miniwdl-lsf.git
```

## Create miniwdl configuration file
In the `miniwdl_setup` directory of this repository, I have included a sample miniwdl configuration file. You will need to place this file in your HOME directory: `~/.config/minidwdl.cfg`

The easiest solution is to clone this repository and move the file into it's final location.

``` 
git clone https://github.com/RhettRautsaw/StJude_PacBio-WDL-tutorial

mkdir -p ~/.config

mv StJude_PacBio-WDL-tutorial/miniwdl_setup/miniwdl.cfg ~/.config/miniwdl.cfg
```

## Test miniwdl installation
To test the miniwdl installation and configuration file, I've also included a small WDL workflow that will scatter 10 jobs onto your HPC, call a base docker/singularity container, and generate 10 "hello_*.txt" files. We will tell miniwdl to run and place the results in `~/WHALE_POD_TEST`
```
miniwdl run StJude_PacBio-WDL-tutorial/miniwdl_setup/whale_pod.wdl --dir ~/WHALE_POD_TEST
```

If you run this command a second time, it should complete much faster as it will locate the cached result from the previous successful run. 


# STAY TUNED FOR MORE
