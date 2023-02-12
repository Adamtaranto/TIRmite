# Tirmite updates

## Update for github features
Add secrets for pypi / conda / gpg key etc
Add issues for dev items
Link issues into milestones
Add codespace settings - Dockerfile etc

## Core functionality
Use bedtools fasta index
Update stderr/stdout reports
Update check for existing files methods
Add progress reporting
Pair with kd tree
## Package structure
Clean up __init__.py files
Split functions into modules
Add if name == main options

## Entry points
tirmite extract - extract TIRs from list of elements
timite build - cluster, align TIRs and build HMM
tirmite find - query genome for HMM hits
tirmite merge - merge two sets of hits
tirmite compare - report unique hits in each set
tirmite pair - take list of hits (HMM or BLAST) and perform pairing
tirmite report - output gffs for elements and TIRs, output fasta, report length distribution, stats on found elements
tirmite classify - cluster complete elements, search for known transposase domains, update gff with cluster labels

// Add two-part entry points https://www.youtube.com/watch?v=0W0k6zP_Lto

## Workflow
Extract TIRs from starting set of elements
Cluster on identity
Align
Build HMMs
Search genome for HMM hits
Search genome for BLAST hits
filter on hit qual / length
Merge overlapping hits in same orientation - bedtools
opt Pair and retain paired hits
extract new TIRs
opt manual trim
Rebuild HMMs
Search for HMM hits (opt with diff model subsets)
opt merge hits
Pair hits
Report annotations
extract complete elements
Cluster complete elements
search for TPase domains

## Tests
Set up for pytest
Add action for pytest
Add test data

## Docs
Add docstring for each function
Add workflow image
Add github pages branch 
// https://docs.github.com/en/pages/getting-started-with-github-pages/creating-a-github-pages-site
// https://docs.github.com/en/pages/setting-up-a-github-pages-site-with-jekyll/adding-a-theme-to-your-github-pages-site-using-jekyll
// https://docs.github.com/en/pages/getting-started-with-github-pages/configuring-a-publishing-source-for-your-github-pages-site

## Publishing
Add OCI/Docker file
Add action for publish to: ghpr, pypi, conda
Publish on new release tag
Add nextflow file