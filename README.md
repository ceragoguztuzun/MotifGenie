# MotifGenie
MotifGenie is a Python command line tool for searching transcription factor binding motif sequences in merged ChIP-Seq binding sites. First, MotifGenie analyzes multiple ChIP-Seq samples simultaneously to identify shared binding regions for a locus of interest, then, searches the binding motif sequence of a given transcription factor in these shared binding regions.

*Authors:* Çerağ Oğuztüzün, Pelin Yaşar, Kerim Yavuz, Mesut Muyan, Tolga Can

![mg1](https://user-images.githubusercontent.com/38559757/117874575-cdfc1900-b2a9-11eb-8677-bf2997f75091.jpg)
*Artwork by:* [brenazan](https://www.instagram.com/brenazan/)


## Workflow  
![MotifGenieWorkflow](https://user-images.githubusercontent.com/38559757/114838020-f85ad380-9ddc-11eb-9643-86b475f34198.png)

MotifGenie is composed of two main modules: 1) finding common binding regions in multiple ChIP-Seq samples and 2) searching for binding sequences using a binding profile of a transcription factor. Peaks identified by MACS 2.0 for 11,286 samples in the Cistrome database are used by the first module. Given a cell line, a TF, a genomic locus, and an occurrence percentage threshold t, the first module initially identifies the subset of samples for the given cell line and the TF. Then, it scans each bed file in that subset for the given genomic locus. The second module uses the JASPAR TF binding profile and searches the genomic sequences in common binding regions. Sequences from the human reference genome, hg38, as provided by the UCSC Genome Browser are used. The top 15 highest scoring sequences in binding regions are given as the output of MotifGenie. A sequence logo of the top-scoring binding sequences is also generated if the corresponding command-line option is used.

## Requirements 
beautifulsoup4-4.9.3
`$ pip3 install bs3` 

logomaker-0.8 
`$ pip3 install logomaker`
## Installation  
Download the mg.py script in a working directory.
## Input Files 
MotifGenie requires a pfm file from the JASPAR database as input, hence the mg.py script must be in the same directory as the pfm file. Pfm file contains a Position Frequency Matrix of the motif to be searched.

An example pfm file of the ESR1 transcription factor from the JASPAR database, where MA0112.1 refers to the ID:
```
>MA0112.1 ESR1
  1.00   1.00   7.00   2.00   0.00   0.00   0.00   6.00   1.00   2.00   3.00   1.00   1.00   5.00   0.00   0.00   2.00   3.00
  5.00   5.00   1.00   0.00   0.00   0.00   7.00   0.00   7.00   5.00   2.00   1.00   0.00   1.00   8.00   9.00   4.00   4.00
  1.00   1.00   1.00   7.00   9.00   0.00   2.00   2.00   1.00   1.00   4.00   1.00   8.00   3.00   0.00   0.00   0.00   1.00
  2.00   2.00   0.00   0.00   0.00   9.00   0.00   1.00   0.00   1.00   0.00   6.00   0.00   0.00   1.00   0.00   3.00   1.00
```
## How To Run
### Parameters
MotifGenie is used by a set of required and optional arguments through the command line interface. MotifGenie does search on the hg38 genome. The cell line information to be inputted must be present in CellLine_list.txt found in the repository. Also, the name of the Transcription Factor pfm file to be inputted must be present in the TF_list.txt file which is found in the repository. If the TF name is not in the list, peaks will not be found.

The **required** arguments include:
- `-jf JASPARFILE` or `--jasparFile JASPARFILE`: Pfm file name from JASPAR. (for example MA0112.1.pfm)
- `-cl CELLLINE` or `--cellLine CELLLINE`: Cell line input. (for example MCF-7)
- `-locus LOCUS` or `--locus LOCUS`: Written as: chr5:1,234,234-1,245,345 or chr5:1234234-1245345
- `-pit PEAKINTERVALTHRESHOLD` or `--peakIntervalThreshold PEAKINTERVALTHRESHOLD`: Threshold of DNA binding sequence to be used in peak locations to find peak sequence intervals. This means filtering the peak sequences where binding percentage passes the inputted threshold value.

The **optional** arguments include:
- `-h` or `--help`: displays the documentation
- `-gene GENENAME` or `--geneName GENENAME`: Gene name to search using UCSC Genome Browser.
- `-ust UPPERSCORETHRESHOLD` or `--upperScoreThreshold UPPERSCORETHRESHOLD`: Upper score threshold for displaying output.
- `-lst LOWERSCORETHRESHOLD` or `--lowerScoreThreshold LOWERSCORETHRESHOLD`: Lower score threshold for displaying output.
- `-n TOPN` or `--topN TOPN`: Number of motifs to be displayed that have the highest scores. (15 by default)
- `-cms` or `--coreMotifSearchFlag`: Runs the search with Core Motifs.
- `-trc THRESHOLDRATIOFORCONSENSUS` or `--thresholdRatioForConsensus THRESHOLDRATIOFORCONSENSUS`: Positions of consensus get trimmed if the frequency ratio is greater than the inputted threshold. Used with coreMotifSearchFlag, when searching using Core Motifs. (0.66 by default)
- `-sl` or `--seqLogoFlag`: Displays and saves the Sequence Logo after motif search when the argument is used.

The command line interface structure: 
```
mg.py [-h] -jf JASPARFILE -cl CELLLINE -locus LOCUS -pit
             PEAKINTERVALTHRESHOLD [-ust UPPERSCORETHRESHOLD]
             [-lst LOWERSCORETHRESHOLD] [-gene GENENAME] [-n TOPN]
             [-trc THRESHOLDRATIOFORCONSENSUS] [-cms] [-sl]
```
### Sample Search Commands
- This searches the motif given through the pfm file name, on the given locus using the cell line and the peak interval threshold (which are the must-have arguments).

`python mg.py -jf MA0112.1.pfm -cl MCF-7 -locus chr5:139638349-139693882 -pit 20`
- This search is similar to the above one but does a core motif search by inputting the threshold ratio for consensus (-trc) value.

`python mg.py -jf MA0112.1.pfm -cl MCF-7 -locus chr5:139638349-139693882 -pit 20 -trc 0.7 -cms`
- This does search on the overlapping sequences of CXXC5 gene's sequences obtained from the UCSC database, and the peak sequences. The search also uses core motif search and specifies the number of outputs (-n) to be 35 motifs.

`python mg.py -jf MA0112.1.pfm -cl MCF-7 -locus chr5:139638349-139693882 -pit 20 -n 35 -gene CXXC5 -trc 0.7 -cms`
- This search is similar to the above one, but specifies the upper score threshold (-ust) and lower score threshold (-lst) to filter the motifs to be outputted. The sequence logo flag (-sl) can be included to each of these commands for viewing the Sequence Logo of the outputted set of motifs.

`python mg.py -jf MA0112.1.pfm -cl MCF-7 -locus chr5:139638349-139693882 -pit 20 -n 30 -gene CXXC5 -ust 8.3 -lst 6.9`
## Output Description 
MotifGenie outputs the search results as a table in the command line interface. If the Sequence Logo flag is used, the generated Sequence Logo will be saved in the directory of the script as a PNG file.

A sample output:
```
SEQUENCE                      SCORE               POSITION                      PARALLELITY     PEAK
-----------------------------------------------------------------------------------------------------
GTGTGCCCTCTACTGGCAG           9.342963            127734105 - 127734123         -               peak1
TTGCGCCACCTGAAGGAGA           8.430384            127737843 - 127737861         -               peak3
CCGCCCCCTCTCTGGGCAG           6.603304            127734156 - 127734174         -               peak1
GCCGCGGGAGGTGGCGGTG           6.159184            127737898 - 127737916         +               peak3
GGGTCAGGTGGGGGCAGGA           6.014313            127734051 - 127734069         +               peak1
CGCCCTTAAGAAGCCGCGG           5.780073            127737910 - 127737928         +               peak3
GCGCCACCTGAAGGAGAAG           5.735097            127737841 - 127737859         +               peak3
GGGGCAGGAGCAGGAGCGT           5.678570            127734041 - 127734059         +               peak1
TCCAGCCACCTCCTTGTTA           5.543466            127733968 - 127733986         -               peak1
TCCACCCTAGCCGGCCGCC           5.353720            127736273 - 127736291         +               peak2
CCACCTGAAGGAGAAGGCG           5.330792            127737838 - 127737856         +               peak3
TTTCCAGCGGGGGAAGGAC           5.273824            127734009 - 127734027         +               peak1
TTGGCTGCAGAAGGTCCGA           5.217327            127733939 - 127733957         +               peak1
CTGCCTCTCGCTGGAATTA           4.891351            127736242 - 127736260         +               peak2
TGGGCGCTAGCGGCTGCGT           4.650896            127737810 - 127737828         +               peak3
```

- `HEADER`: the summary of the search parameters such as, whether a specific gene name was specified, core motif search was used and the inputted consensus trimming threshold.
- `SEQUENCE`: the sequence of motifs found sorted with respect to the scores.
- `SCORE`: the log-based score given to the motif which is correlated with the similarity of the inputted motif that was searched on the reference sequence.
- `POSITION`: the start and end locations where the motif was found on the reference genome.
- `PARALLELITY`: whether the motif was found in a parallel (+) or antiparallel (-) search.
- `PEAK`: the number of peak sequence where the motif was found on.
- `ENTRY`: the entry number referring to the UCSC database, in which the sequence for the inputted gene was found. This column is only presented for gene-specific search.

Sample Sequence Logo output of the above result
![MA01391](https://user-images.githubusercontent.com/38559757/109424770-33a06f00-79f6-11eb-83ec-fdbf9fc6b904.png)

