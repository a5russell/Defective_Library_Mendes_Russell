# Defective_Library_Mendes_Russell
Repository for code on our analysis of defective influenza species using a length-based library analysis
Only files smaller than 10mb were pushed to this repository. 
Larger files, such as sequencing, will need to be downloaded and reconstructed in the apppropriate folders in order to re-run this analysis.

Code was run with the following accessible from the command line (in PATH)
Samtools 1.9
HTseq 0.11.3
STAR 2.7.1a
BLASTn 2.9.0+
Trimmomatic 0.39
Python 3.7
R 4.0.2

and the following python packages
numpy 1.19.5
matplotlib 3.3.3
seaborn 0.11.0
pandas 1.2.0
scipy 1.6.0
statsmodels 0.12.1

and the following R package
DESeq2 1.28.1

Scripts files contains utility scripts.

Analysis is split across four jupyter notebooks and two R markdown files.

Processing_and_analysis_of_barcoded_libraries.ipynb describes the generation and analysis of barcoded libraries.
Interferon_beta_natural_diversity.ipynb describes the analysis of naturally-occurring deletions.
qPCR_analysis.ipynb describes the processing and graphing of all qPCR data.
Flow_cytometry.ipynb describes the processing and graphing of all flow cytometry data after initial debris gate set in FlowJo and data exported to csv. 

expressionAnalysis.Rmd. Combined with Interferon_beta_natural_diversity.ipynb analyzes mRNASeq data validating interferon sort.
barcode_DESeq.Rmd. Combined with Processing_and_analysis_of_barcoded_libraries.ipynb short script to analyze barcoded junctions using DESeq2.

The following folders contain data or analyses.
Database- contains WSN and human genomes necessary for analysis.
Results- Results of analyses.
qPCR- Data for qPCR analyses.
flowCytometry - Data for flow cytometry analyses
Sequencing- Folders that would contain sequencing data to reconstruct analysis.
