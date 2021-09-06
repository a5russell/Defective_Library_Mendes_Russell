# Defective_Library_Mendes_Russell

The code contained within this repository describe the analysis of data, and generation of figures, in our attempts to understand length, segment, and replication-dependent phenomena of defective influenza species.

Analysis are split across several different jupyter notebooks, with descriptions below, and figures inline.

All intermediate files less than 10mb were pushed to this repository to aid in understanding our analysis.
Larger files, such as sequencing, will need to be downloaded and reconstructed in the apppropriate folders in order to re-run this analysis.
Specifically, place NGS files within appropriate directories within Sequencing to rerun all sequencing analyses.

## Directories

The following directories, and their purpose, exist within this repository.

- <b>Database</b>       Repository for human and influenza genomic sequences. A/WSN/1933 BLAST database and STAR indicies provided. To regenerate the human genome,                         Gchr38 must be concatenated to A/WSN/1933 and star indices generated from the combined genome.
- <b>Results</b>        Final datafiles for analyses after processing. Most provided as supplemental files within manuscript.
- <b>Scripts</b>        Short scripts written for this analysis. Seperated from jupyter notebooks for readability and portability.
- <b>Sequencing</b>     Folder that would contain NGS samples to regenerate this pipeline. Folder achitecture essential to run of pipeline.
- <b>flowCytometry</b>  Flow cytometry raw data
- <b>qPCR</b>           qPCR raw data
  
  

## Dependencies

Code within this repository was run with the following tools, and versions, installed and available from PATH. Specific websites and documentation are provided where available. 

- <b>Python</b>      run with version 3.7. Available from https://www.python.org/downloads/release/python-370/
- <b>Trimmomatic</b> run with version 0.39. Available from http://www.usadellab.org/cms/?page=trimmomatic
- <b>HTseq</b>       run with version 0.11.3. Available from https://htseq.readthedocs.io/en/master/
- <b>STAR</b>        run with version 2.7.1.a. Available from https://github.com/alexdobin/STAR
- <b>Samtools</b>    run with version 1.9. Available from http://www.htslib.org/
- <b>FastQC</b>      run with version 0.11.8. Available from https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
- <b>BLASTn</b>      run with version 2.9.0+. Available from https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download
- <b>R</b>           run with version 4.0.2. Available from https://cran.r-project.org/bin/windows/base/

The following python packages and versions were used. All were installed using conda. (https://docs.conda.io/en/latest/)
- <b>numpy</b>       run with version 1.19.5. (https://numpy.org/)
- <b>matplotlib</b>  run with version 3.3.3. (https://matplotlib.org/)
- <b>seaborn</b>     run with version 0.11.0. (https://seaborn.pydata.org/)
- <b>pandas</b>      run with version 1.2.0.(https://pandas.pydata.org/)
- <b>scipy</b>       run with version 1.6.0. (https://www.scipy.org/)
- <b>statsmodels</b> run with version 0.12.1. (https://www.statsmodels.org/stable/index.html)

The following R packages and versions were used. All were installed using CRAN (https://cran.r-project.org/)
- <b>DESeq2</b>      run with version 1.28.1



