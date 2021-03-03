# TranSuite

A software suite for accurate identification, annotation, translation, and feature characterization of annotate transcripts.
Reference: 
* Juan C Entizne, Wenbin Guo, Cristiane P.G. Calixto, Mark Spensley, Nikoleta Tzioutziou, Runxuan Zhang, and John W. S. Brown; [*TranSuite: a software suite for accurate translation and characterization of transcripts*](https://www.biorxiv.org/content/10.1101/2020.12.15.422989v1) (*Peer review in progress*).


----------------------------
# Table of Contents
----------------------------

   * [Overview](#overview)
   * [Installation](#installation)
   * [Modules](#modules)
   * [Input files](#input-files)
   * [FindLORF](#findlorf)
      * [Command and options](#command-and-options)
      * [Output files](#output-files)
   * [TransFix](#transfix)
      * [Command and options](#command-and-options-1)
      * [Output files](#output-files-1)
   * [TransFeat](#transfeat)
      * [Command and options](#command-and-options-2)
      * [Output files](#output-files-2)
   * [Auto](#auto)
      * [Command and options](#command-and-options-3)
      * [Output files](#output-files-3)
   * [Future work](#future-work)
   * [License](#license)
   * [Contact](#contact)


----------------------------
# Overview
----------------------------

TranSuite is a software for identifying coding sequences of transcripts, selecting translation start sites at gene-level, generating accurate translations of transcript isoforms, and identifying and characterizing multiple coding related features, such as: coding potential, similar-translation features, and multiple NMD-related signals. TranSuite consists of three independent modules, FindLORF, TransFix and TransFeat (Fig. 1). Each module can be run independently or as a [pipeline with a single command](#auto).

- FindLORFS - finds and annotates the longest ORF of each transcript
- TransFix - “fixes” the same translation start codon AUG in all the transcripts in a gene and re-annotates the resulting ORF of the transcripts (Fig. 2)
- TransFeat - identifies structural features of transcripts, coding potential and NMD signals (Fig. 3)
- Auto - executes the whole pipeline (FindLORFS, TransFix, and TransFeat) in tandem (Fig. 1)


![Fig1_TS_pipeline.jpg](https://drive.google.com/uc?id=1V6AEorw85_VGCne6Na9uSaWUWeg-a0Y8)

**Fig 1.**  TranSuite pipeline. The TranSuite pipeline takes A) the transcript annotations to be analysed (in GTF format) and B) the transcript sequence (FASTA file) as inputs. 1) FindLORF - For each transcript the longest ORF is identified and genomic co-ordinates generated. 2) TransFix – all transcripts of a gene are translated with the fixed translation start AUG. 3) TransFeat – transcript features are derived and transcripts classified. Information on genes and transcripts is contained in a TranSuite-generated report and tables.


![Fig3_TransFix_example.jpg](https://drive.google.com/uc?id=1ZP22Olbg_Zr4fLPkhWRKA-YlphDMQvSh)

**Fig 2.** TransFix correctly translates and annotates ORFs in transcript variants. A) Schematic of gene with a fully spliced (FS) transcript that encodes the full-length protein and an alternatively spliced (AS) transcript where the AS changes the frame of translation and introduces a PTC (TransFix). B) Common mis-annotation of the AS transcript due to selection of the longest ORF by translation programs. This selects an AUG downstream of the PTC (ignoring the authentic translation start site) producing a long downstream ORF (ldORF) which represents a C-terminal fragment of the protein of the gene (top). Where an AUG in one of the other two reading frames occurs in frame and upstream of the ldORF, the ORF is extended (grey boxes). The ORF contains a C-terminal fragment of the protein of the gene and an unrelated N-terminal fragment. White boxes – protein-coding exons; black boxes – UTRs; grey boxes – coding exons with sequence from different frame; arrow – position of authentic translation start AUG.


![Fig2_TransFeat_workflow.jpg](https://drive.google.com/uc?id=1_Lasc18vsbYccNHJ8pG6xXCf7-DMVfxM)

**Fig 3.** Transcript feature characterisation by TransFeat. A) translation data identifies transcripts with non-coding features (no CDS - no AUG; CDS < 30 amino acids or short ORF where CDS is between 30-100 amino acids (default values); Genes with only no CDS transcripts are classed as Non-coding genes; B) All other genes and their associated transcripts are protein-coding genes. Transcripts containing a PTC (PTC+) are unproductive and all other transcripts are protein coding isoforms. Transcripts of some genes code for identical proteins if AS occurs only in the 5’ and/or 3’UTRs. Transcripts with NAGNAG alternative splicing code for proteins which differ by one amino acid. Other protein-coding transcripts encode protein variants. C) PTC+ transcripts are further characterised by different features: downstream splice junction, long 3’UTR and overlapping uORF (NMD signals), in frame and out of frame uORFs, long downstream ORFs. The dotted line from in/out of frame uORFs reflects the potential of some uORFs to trigger NMD.


----------------------------
# Installation
----------------------------

TranSuite has been developed in Python 3.6 

TranSuite requires the following packages:
- [BioPython v1.78](https://anaconda.org/anaconda/biopython)


Packages installations commands:
```
conda install -c anaconda biopython
```

TranSuite is ready to use. The compressed file can be directly downloaded from the [GitHub repository](https://github.com/anonconda). Once uncompressed, TranSuite can be used directly from the command line by specifying the path to the main executable `transuite.py`

Additionally, in the future it will also be possible to install TranSuite through popular Python installation managers PyPI and Anaconda:


```


----------------------------
# Modules
----------------------------

TranSuite works with a module / module-options structure:

```
transuite.py module options
```
where the modules are:

- **FindLORF**    : Finds and annotates the longest ORF of each transcript.
- **TransFix**       : Fix the same translation start codon AUG in all the transcripts in a gene and re-annotates the resulting ORF of the transcripts.
- **TransFeat**       : Identify structural features of transcripts, coding potential and NMD signals.
- **Auto**        : Execute the whole pipeline (FindLORFS, TransFix, and TransFeat) in tandem.

Each module can be executed as follow:

```
python /path/to/transuite.py module options
```

For example, to observe the help documentation of Auto module:

```
python /path/to/transuite.py Auto --help
```


----------------------------
# Input files
----------------------------

All of TranSuite modules (FindLORF, TransFix and TransFeat) use the same input files format:

- The transcriptome annotation to analyze, in [GTF format](https://www.ensembl.org/info/website/upload/gff.html)
- The transcripts nucleotide sequences, in [FASTA format](http://bioinformatics.org/annhyb/examples/seq_fasta.html)

Please note that the programs assume that nucleotide sequences in the FASTA file represent the exonic region of the transcript. Please beware that errors will happen if the user provide a FASTA file of the coding-region (CDS) sequences instead. 


Additional notes:

- FindLORF will parse only the "exon" [feature](https://www.ensembl.org/info/website/upload/gff.html) information from the GTF file
- TransFix and TransFeat will extract "CDS" [feature](https://www.ensembl.org/info/website/upload/gff.html) information from the GTF file
- Transcripts without annotated "CDS" features will remain unprocessed by TransFix, and they will be identified as "No ORF" by TransFeat
- When executing the whole pipeline (Auto module), TranSuite will automatically forward the appropiate GTF file to each module. That is: FindLORF output GTF will be use as input by TransFix, and TransFix output GTF will be TransFeat input


----------------------------
**FindLORF**
==============

----------------------------

FindLORF identifies and annotate ORF information in newly transcriptome annotations. Firstly, FindLORF translates each transcript sequence in its three frames of translation according to its annotated strand and stores the relative start and stop codon positions of all the resulting ORFs. Secondly, FindLORF selects the longest ORF for each transcript as its putative CDS region. Finally, FindLORF annotates the CDS using the genomic information contained in the transcriptome annotation to convert the relative ORF start-stop codon positions in the transcript sequence into genomic co-ordinates. The FindLORF module takes as input the transcriptome annotation to be curated (GTF format) and the transcripts exon sequences (FASTA format). See [Input files](#input-files) above.

## Command and options
Command to run FindLORF:
```
python transuite.py FindLORF [options]
```
```
python transuite.py FindLORF --gtf <input-gtf.gtf> --fasta <input-fasta.fa> --outpath </path/for/output-folder> --outname <outname> --cds <30>
```

List of options available:
- **--gtf**: Transcriptome annotation file in GTF format
- **--fasta**: Transcripts nucleotide fasta file
- **--outpath**: Path of the output folder
- **--outname**: Prefix for the output files
- **--cds**: Minimum number of amino-acids an ORF must have to be considered as a potential CDS. Default: 30 AA

Example:
```
python transuite.py FindLORF --gtf ./test_dataset/subset_AtRTD2.gtf --fasta ./test_dataset/subset_AtRTD2_transcripts.fa --outpath ./test_dataset/test_output --outname  --cds 30
```

## Output files
FindLORF automatically generates a subfolder to store the output files:
/**&lt;outpath&gt;**/**&lt;outname&gt;**\_longorf/

FindLORF generates the following output files:
1. *GTF* file with the longest ORF in the transcripts annotated as its CDS
2. *FASTA* files of the transcripts CDS regions (nucleotide, and peptide)
3. Log *CSV* files reporting transcripts that could not be annotated, for example: due to lack of an AUG
4. A *JSON* file containing the transcripts ORF relative coordinates


-------------------
**TransFix**
==============
-------------------

TransFix provides more biologically accurate translations by selecting the authentic translation start site for a gene, “fixing” this location and using it to translate the gene transcripts and annotating the resulting CDS of the translations. We define the authentic translation start site as the site used to produce the full-length protein of the gene. In detail, TransFix firstly extracts the CDS co-ordinates of the transcripts from the transcriptome annotation and groups the transcripts according to their gene of origin. Then, TransFix selects the start codon of the longest annotated CDS in the gene as the representative translation start site and translates all of the transcripts in the gene from the “fixed” translation start site. Finally, TransFix annotates the genomic co-ordinates of the resulting stop codons. In some cases, transcript isoforms do not contain the “fixed” translation start site due to an AS event or an alternative transcription start site. To account for this, TransFix tracks those transcripts that are not translated during the first fix AUG/translation cycle and they are then processed through a second fix AUG/translation cycle to determine and annotate their valid translation start-sites.

### **Command and options** ###
Command to run TransFix:
```
python transuite.py TransFix [options]
```
```
python transuite.py TransFix --gtf <input-gtf.gtf> --fasta <input-fasta.fa> --outpath </path/for/output-folder> --outname <outname> --iter <5>
```

List of options available:
- **--gtf**: Transcriptome annotation file in GTF format
- **--fasta**: Transcripts nucleotide fasta file
- **--outpath**: Path of the output folder
- **--outname**: Prefix for the output files
- **--iter**: Maximum number of 'start-fixing & translation' cycles to identify alternative start sites. Default: 5
- **--chimeric**: Table indicating chimeric genes in the annotation (Optional) 

Example:
```
python transuite.py TransFix --gtf ./test_dataset/test_output/test_run_longorf/test_run_longorf.gtf --fasta ./test_dataset/subset_AtRTD2_transcripts.fa --outpath ./test_dataset/test_output --outname test_run --iter 5
```

## Output files
TransFix automatically generates a subfolder to store the output files:
/**&lt;outpath&gt;**/**&lt;outname&gt;**\_transfix/

TransFix generates the following output files:
1. *GTF* file with the fixed CDS coordinates at the gene-level
2. *FASTA* files of the transcripts CDS regions (nucleotide, and peptide)
3. Multiple log *CSV* files: a) log files reporting transcripts that could not be annotated, for example for lack of an annotated CDS; and b) logfile tracking the *fixing* cycle at which the AUG was annotated


-------------------
**TransFeat**
==============
-------------------

TransFeat extracts and processes the transcripts CDS information contained in transcriptome annotations to infer multiple characteristics of the genes, transcripts and their coding potential (Fig. 2), and it to reports theis information in an easily accessible format.

### **Command and options** ###
Command to run TransFeat:
```
python transuite.py TransFeat [options]
```
```
python transuite.py TransFeat --gtf <input-gtf.gtf> --fasta <input-fasta.fa> --outpath </path/for/output-folder> --outname <outname> --pep <30> --ptc <70>
```
List of options available:
- **--gtf**: Transcriptome annotation file in GTF format
- **--fasta**: Transcripts nucleotide fasta file
- **--outpath**: Path of the output folder
- **--outname**: Prefix for the output files
- **--pep**: Minimum number of amino-acids a translation must have to be consider a peptide. Default: 100 AA
- **--ptc**: Minimum CDS length percentage below which a transcript is considered prematurely terminated (PTC). Default: 70%

Example:
```
python transuite.py TransFeat --gtf ./test_dataset/test_output/test_run_transfix/test_run_transfix.gtf --fasta ./test_dataset/subset_AtRTD2_transcripts.fa --outpath ./test_dataset/test_output --outname test_run --pep 30 --ptc 70
```

**Note:** When running the analysis on the above *test_dataset* you will get a *WARNING* message regarding a number of features not present in the TransFeat table. This is expected given the small number of transcripts in the test dataset.

### **Output files** ###
TransFeat automatically generates a subfolder to store the output files:
/**&lt;outpath&gt;**/**&lt;outname&gt;**\_longorf/

TransFeat generates the following output files:
1. A main *CSV* table reporting the transcripts coding features
2. *FASTA* files of: (**1**) transcripts with CDS (both protein-coding and unproductive), (**2**) transcripts classified by coding-potentiality (protein-coding transcripts, non-coding genes), (**3**) transcripts alternative ORFs (uORF, ldORF) 
3. Multiple *CSV* tables reporting number of transcripts and transcripts-features subdivided by gene categories
4. Multiple *CSV* tables reporting co-ordinates and sequences of transcripts subdivided by feature categories (ldORF, NMD)
5. A *JSON* file containing the transcripts ORF relative coordinates


-------------------
**Auto**
==============
-------------------

This module performs FindLORF, TransFix, and TransFeat analysis in tandem with a single command. 

### **Command and options** ###
Command to run Auto:
```
python transuite.py Auto [options]
```
```
python transuite.py Auto --gtf <input-gtf.gtf> --fasta <input-fasta.fa> --outpath </path/for/output-folder> --outname <outname> --cds <30> --iter <5> --pep <100> --ptc <70>
```
List of options available:
- **--gtf**: Transcriptome annotation file in GTF format
- **--fasta**: Transcripts nucleotide fasta file
- **--outpath**: Path of the output folder
- **--outname**: Prefix for the output files
- **--cds**: Minimum number of amino-acids an ORF must have to be considered as a potential CDS. Default: 30 AA
- **--iter**: Maximum number of 'start-fixing & translation' cycles to identify alternative start sites. Default: 5
- **--pep**: Minimum number of amino-acids a translation must have to be consider a peptide. Default: 100 AA
- **--ptc**: Minimum CDS length percentage below which a transcript is considered prematurely terminated (PTC). Default: 70%
- **--chimeric**: Table indicating chimeric genes in the annotation (Optional)

Example:
```
python transuite.py Auto --gtf ./test_dataset/subset_AtRTD2.gtf --fasta ./test_dataset/subset_AtRTD2_transcripts.fa --outpath ./test_dataset/test_output --outname test_run --cds 30 --iter 5 --pep 100 --ptc 70
```

### **Output files** ###
The Auto module automatically generates all of the modules subfolders and their respective output files:
/**&lt;outpath&gt;**/**&lt;outname&gt;**\_longorf/
/**&lt;outpath&gt;**/**&lt;outname&gt;**\_transfix/
/**&lt;outpath&gt;**/**&lt;outname&gt;**\_transfeat/

The main output files of TranSuite pipeline are:
1. The *GTF* file generate by TransFix
2. The main *CSV* feature table generate by TransFeat
3. The *FASTA* files generate by TransFix and/or as classified by TransFeat (Coding transcripts, Non-coding genes)
4. Any of the multiple log files generated during the analysis


----------------------------
# Future work
----------------------------

TODO:
- Create PyPI / Anaconda installer
- Allow user to modify the minimum length to identify uORF, and Long 3'UTR (Current values: 10 AA and 350 nt respectively)
- Allow user to modify the minimum "nucleotide distance from last Splice-Junction" for the identification of the DSSJ NMD-signal (default 50 nt)
- Tidy up and make available *in-house* scripts to identify overlapping gene IDs in the annotation (This script could be useful to generate tables of potentially chimeric gene IDs)
- Add a more detailed description of TransFeat output table into the README
- Create a Web-Tool interface for TranSuite in [Galaxy](https://usegalaxy.org/)
- Allow the use of alternate codon tables


----------------------------
# License
----------------------------

TranSuite is released under the [MIT license](https://opensource.org/licenses/MIT)


----------------------------
# Contact
----------------------------

For any further enquiries please contact the main developer at <e.entizne@dundee.ac.uk>
