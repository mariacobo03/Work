---
### BEGIN-Of-YAML-Block ###
#
## ######################################################################################
##
##   README_BScBICG2425_exercise04_SURNAME_NAME.md
##
##   A LaTeX-extended MarkDown template for BScBI-CG practical exercise submissions.
##
## ######################################################################################
##
##                 CopyLeft 2024 (CC:BY-NC-SA) --- Josep F Abril
##
##   This file should be considered under the Creative Commons BY-NC-SA License
##   (Attribution-Noncommercial-ShareAlike). The material is provided "AS IS", 
##   mainly for teaching purposes, and is distributed in the hope that it will
##   be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
##   of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
##
## ######################################################################################
#
# The current execise number
thyexercise: 04
#
# title-meta is title string without any LaTeX format to be used as pdftitle, part of emails subject...
title-meta: BScBI-CG Exercise 04
#
# title is the big title for the cover page, fully LaTeX formated to fit into a shortstack command...
title: |
  \newline \textsc{BScBI-CG}
  \newline \textsc{Practicals}
  \newline \textsc{Report}
#subtitle:
#
# runtitle is the running header or footer, used i.e. by fancyheadings...
runtitle: |
  BScBI-CG Exercise 04 Report
#
# author-meta sets the pdfauthor variable...
author-meta: !!str 'Name SURNAME @ BScBI Computational Genomics'
#
# authors to appear on the title page...
author:
- name: Maria Cobo
  myemail: maria.cobo
  mydomain: alum.esci.upf.edu
#
# IMPORTANT: we need to define email as two fields (user/domain)
#            to avoid parsing email by YAML, as in the following example.
#
# - name: Josep F Abril
#   myemail: jabril   #
#   mydomain: ub.edu  # the real complete email address was in this case: jabril@ub.edu
#
# authorshort defines a brief list of authors for headings, i.e.: Abril, JF
authorshort: Cobo, M
#
# template formating variables...
papersize: a4paper
fontsize: 10pt
geometry: margin=1.5cm
toc: true
lof: true
lot: true
colorlinks: true
urlcolor: blue
citecolor: green
linkcolor: red
#
# LaTeX definitions (available for the main document)
further-defs: 
-  \def\CGVC{\href{https://aula.esci.upf.edu/course/view.php?id=8516}{Computational Genomics Virtual Campus at ESCI}}
-  \def\ScerS{\textit{Saccharomyces cerevisiae} (strain S288C)}
-  \def\scerS{\textit{S.~cerevisiae} (strain S288C)}
-  \def\Scer{\textit{Saccharomyces cerevisiae}}
-  \def\scer{\textit{S.~cerevisiae}}
-  \def\sra{\textsc{SRA}}
-  \def\Sra{\textsc{Short Reads Archive}}
-  \def\SRA{\href{https://www.ncbi.nlm.nih.gov/sra}{\sra}}
#
### End-Of-YAML-Block ###
---

```{=comment}
We do not need the comment LaTeX environment to hide a block,
with pandoc >2.x one can define "comment" blocks, which are even
escaped and not passed to LaTeX and pandoc (achieving the same hidden effect).
\begin{comment}
The \texttt{comment} \LaTeX\ environment allow us to include any piece of text, 
code, etc, but it wil not be included into the final PDF report when compiling 
the \texttt{MarkDown} file. You can open/close one of such environments 
at any time if you need them.
\end{comment}
```

# Introduction

From the previous exercise, we are going to check completeness of the
*de novo* assemblies we have generated for \Scer\ and to compare them
against the available reference genome. To speed up those comparisons
we can focus on specific chromosomes, so that we will need to filter
out sequences or alignments to continue with the analyses. After that
we can take a look to the repetitive elements that may be present in
contigs from our best assembly set.


## Objectives

* To compare our assembly against a reference genome, so we can assess
  the performance of the assembler and the protocol.
* To check completeness of our assemblies, so that we can choose the
  best assembly for the downstream analyses, such as masking,
  gene-prediction, and so on.
* To complete the assembly protocol by masking the repeats detected on
  the chosen contigs/scaffolds.
* We introduce some \LaTeX\ examples for citing paper references as footnotes.

## Prerequisites

### Installing required software

As for the previous practical, we must ensure that at least `pandoc`
and `pdflatex` commands are running smoothly over our report files. If
you still need to install the base software, please refer to
`exercise_00` and `exercise_01`, as well as the short tutorials from
the \CGVC. Remind that we assume that you are using a Debian-based
linux distribution, so we will show only the corresponding set of
commands for that distribution.

For this practical you probably may need to install the following packages:

```{.sh}
#################################
# We are assuming that you have already installed from previous exercises
# the following packages:
#   samtools bamtools picard-tools bwa bowtie2 igv

# ncbi-blast+ - next generation suite of BLAST sequence search tools
sudo apt-get install ncbi-blast+

# gnuplot-qt - a portable command-line driven interactive data and function plotting utility
#              It is required for mummerplot later on...
sudo apt-get install gnuplot-qt
sudo ln -vfs /usr/bin/gnuplot-qt /usr/bin/gnuplot
# just in case mummerplot returns error: Inappropriate ioctl for device
sudo ln -vfs /usr/bin/gnuplot-qt /etc/alternatives/gnuplot

# graphicsmagick - collection of image processing tools (replacement for imagemagick)
sudo apt-get install graphicsmagick

# mummer - Efficient sequence alignment of full genomes
sudo apt-get install mummer --reinstall
#
# IMPORTANT: Some fixes are needed prior to run gnuplot within mummerplot
#            (install mummer package first)
#
#   + just in case mummerplot returns error: Can't use 'defined(%hash)'
sudo sed -i 's/defined (%/(%/
            ' /usr/bin/mummerplot
#
#   + gnuplot fails because some instruction was implemented in later versions
sudo sed -i 's/^\(.*set mouse.*\)$/#\1/
            ' /usr/bin/mummerplot
#
#   + mummerplot cannot find gnuplot:
sudo sed -i 's/system (\"gnuplot --version\")/system (\"\/usr\/bin\/gnuplot --version\")/
            ' /usr/bin/mummerplot
sudo sed -i 's/my \$cmd = \"gnuplot\";/my \$cmd = \"\/usr\/bin\/gnuplot\";/
            ' /usr/bin/mummerplot
#
#   + just in case mummerplot returns error: Inappropriate ioctl for device
sudo ln -vfs /usr/bin/gnuplot-qt /etc/alternatives/gnuplot;
#
# Further NOTES:
#   You probably need to run all those fixes if your system has already gnuplot version > 4.
#   MacOS users lacking sed can try with "perl -i -pe" instead of "sed -i".

# hmmer2 - profile hidden Markov models for protein sequence analysis
sudo apt-get install hmmer

#
# Note: Before installing BUSCO from source it is advisable to install bbmap and metaeuk
#
# bbmap - BBTools genomic aligner and other tools for short sequences
cd $BIN
wget https://sourceforge.net/projects/bbmap/files/BBMap_39.01.tar.gz/download \
     -O BBMap_39.01.tar.gz
tar -zxvf BBMap_39.01.tar.gz bbmap/
#
# metaeuk - toolkit for large-scale gene discovery
#           and annotation in eukaryotic metagenomic contigs. 
cd $BIN
git clone https://github.com/soedinglab/metaeuk.git
mv metaeuk MetaEuk
cd MetaEuk/
mkdir build
cd build/
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=. ..
make -j
make install
cd $BIN
ln -vs ./MetaEuk/build/bin/metaeuk ./metaeuk
#
# BUSCO - estimating the completeness and redundancy of processed genomic data
#         based on universal single-copy orthologs.
#
# Note: you can also download and install BUSCO using conda,
#       see further instructions at:  https://anaconda.org/bioconda/busco
#
#       If using conda, you can try:
#             conda create -n busco5 -c conda-forge -c bioconda busco=5.4.3
#       Then, run it after activating the environment:
#             conda activate busco5
#       And remember to set up all project vars again as you enter into a new shell.
#
cd $BIN/
git clone https://gitlab.com/ezlab/busco.git
cd busco/
sudo python3 setup.py install
cd $WDR
# 
# wget https://gitlab.com/ezlab/busco/-/archive/5.4.3/busco-5.4.3.tar.gz \
#      -O busco-5.4.3.tar.gz
# tar -zxvf busco-5.4.3.tar.gz
# cd busco-5.4.3
# python3 setup.py install --user
# cd $WDR
# ln -vs ./busco-5.4.3/bin/busco $BIN/busco
# ln -vs ./busco-5.4.3/bin/busco ~/busco --> this worked for me 

# 
# repeatmasker - finds repeat families from biological sequences
#
### IMPORTANT NOTE ### -----------------------------------------------------------------------
#
# As RepeatMasker may have some installation issues due to system and/or conda libraries,
# it is recommended to use the docker containers as explained on the short tutorials
# from the virtual campus. Anyway, it should be easy to install with the following instructions:
#
# Repeatmasker may require some extra dependencies
# and you only need those provided by the commands below.
# For further info you can visit: http://www.repeatmasker.org/RepeatMasker/
#
sudo perl -MCPAN -e'install Text::Soundex'
sudo apt install python3-h5py
#
#   we need to download the binaries for tandem repeats finder from:
#              http://tandem.bu.edu/trf/trf409.linux64.download.html
cd $BIN
git clone https://github.com/Benson-Genomics-Lab/TRF.git
cd TRF
mkdir build
cd build
../configure
make
cd $BIN
ln -vs ./TRF/build/src/trf ./
#
# Now it's time to install RepeatMasker (take care, it's a large file, >800Mbp).
wget http://www.repeatmasker.org/RepeatMasker/RepeatMasker-4.1.3-p1.tar.gz \
     -O $BIN/RepeatMasker-4.1.3-p1.tar.gz
cd $BIN
tar -zxvf RepeatMasker-4.1.3-p1.tar.gz
mv RepeatMasker RepeatMasker-4.1.3
cd RepeatMasker-4.1.3c
sudo perl ./configure
# Important: provide here paths only to trf and to hmmer,
#            and set hmmer as default search engine.
cd $BIN
ln -vs ./RepeatMasker-4.1.3/RepeatMasker .
cd $WDR

### NOTE ###
## You can install docker server and download the repeatmasker container
## following the instructions from the Virtual Campus if you experience problems
## when you try to compile or install this tool.
```


#### Using `conda/mamba` environments:

As we saw in the previous exercises, another way to install the
software required to complete the exercises is to use `conda`
environments. You can install `conda` following the instructions from
[this link](https://conda.io/projects/conda/en/latest/user-guide/install/index.html);
you can also use `mamba` instead, which is a compact and faster
implementation of `conda`, from the instructions at [this
link](https://github.com/conda-forge/miniforge#install). Once you have
one of those environment managers installed, you can follow the
commands in the next code block to create the `BScBI-CG2425_exercises`
environment and activate it. __You probably have the conda environment
created from the previous exercise, then you can jump to the next
block of code.__

```{.sh}
#
# ***Important***: ensure that you run the create command
#                  outside any other environment (even the `base` one),
#                  for a fresh install of the proper dependencies.
#
# If you have conda instead of mamba already installed on your system
# you can just replace 'mamba' by 'conda' on the commands below:
mamba env create --file environment.yml

# Now you can run the tools installed on that environment by activating it:
mamba activate BScBI-CG2425_exercises

# Remember that each time you deactivate a conda environment
# all shell variables defined inside will be lost
# (unless they were exported before activating the conda environment).
# Anyway, you can reload project vars with:
source projectvars.sh

# To return to the initial terminal state, you must deactivate the environment:
mamba deactivate
```

__IMPORTANT:__ For this exercise we only need to update our
environment, in order to include the tools introduced to complete
current the protocol (basically adding `emboss` suite to the current
environment). The `environment.yml` file included in the exercise
tarball is the same as that of `exercise_00`, including an extra
dependency line.

```{.sh}
#
# ***Important***: ensure that you run the update command
#                  outside any mamba/conda environment too.
#
# Again, if you have conda instead of mamba already installed on your system
# you can just replace 'mamba' by 'conda' on the commands below:
mamba env update --file environment.yml

# Now you can run the tools installed on that environment by activating it:
mamba activate BScBI-CG2425_exercises

# Remember that each time you deactivate a conda environment
# all shell variables defined inside will be lost
# (unless they were exported before activating the conda environment).
# Anyway, you can reload project vars with:
source projectvars.sh

# To return to the initial terminal state, you must deactivate the environment:
mamba deactivate
```

You can review the contents of the environment YAML file at the
Appendices (see section \ref{prg:environmentYML} on page
\pageref{prg:environmentYML}). Finally, you may need to
perform some fixes on the mummerplot script that was installed
in your conda environment:

```{.sh}
#
# IMPORTANT: Some fixes are needed prior to run gnuplot within mummerplot
#            (install mummer package first)
#
# Remind that this folder can change form one system to another:
# MUMMER=$HOME/.conda/envs/BScBI-CG2425_exercises/bin/mummerplot;
MUMMER=$HOME/miniconda3/envs/BScBI-CG2425_exercises/bin/mummerplot;
#
#   + just in case mummerplot returns error: Can't use 'defined(%hash)'
sed -i 's/defined (%/(%/' $MUMMER;
#
#   + gnuplot fails because some instruction was implemented in later versions
sed -i 's/^\(.*set mouse.*\)$/#\1/' $MUMMER;
#
#   + mummerplot cannot find gnuplot:
sed -i 's/system (\"gnuplot --version\")/system (\"\/usr\/bin\/gnuplot --version\")/
       ' $MUMMER;
sed -i 's/my \$cmd = \"gnuplot\";/my \$cmd = \"\/usr\/bin\/gnuplot\";/
       ' $MUMMER;
```


### Initializing the main report files

As in the previous exercises, remember to download first the exercise
tarball from the \CGVC, unpack this file, modify the files accordingly
to the user within the exercise folder, and set it as the current
working directory for the rest of the exercise...

```{.sh}
# You probably have already done this step.
tar -zxvf BScBI_CG2425_exercise_04.tgz
cd exercise_04

# Rename report file including your "NAME" and "SURNAME"
mv -v README_BScBICG2425_exercise04_SURNAME_NAME.md \
      README_BScBICG2425_exercise04_yourSurname_yourName.md

# Open exercise files using your text editor of choice
# (for instance vim, emacs, gedit, sublime, atom, ...);
# fix "NAME" and "SURNAME" placeholders on them
# and save those changes before continuing.
emacs  projectvars.sh \
       README_BScBICG2425_exercise04_yourSurname_yourName.md &

# Let's start with some initialization.
source projectvars.sh
echo $WDR

# Once you have run the commands that are already in the initial
# MarkDown document, you are probably ready to run this:
runpandoc
```

Let's start with the analyses, and may the shell be with you...


## Datasets

On the previous exercise we have started to work on genomic datasets
for the budding yeast
\href{https://www.ncbi.nlm.nih.gov/genome/?term=txid559292[Organism:noexp]}{\Scer},
one of the major model organisms for understanding cellular and
molecular processes in eukaryotes. We also told on that exercise to
keep the files for re-using them on the current practical, so that we
assume you already have some of the initial files. The \scerS\
assembly R64 reference genome version\footnote{Engel et al. "The
Reference Genome Sequence of \textit{Saccharomyces cerevisiae}: Then
and Now".\newline\hspace*{1cm} \textit{G3 (Bethesda)},
g3.113.008995v1, 2013
(\href{https://www.ncbi.nlm.nih.gov/pubmed/?term=24374639}{PMID:24374639})}
has 16 chromosome sequences, plus the mitochondrion genome, totaling
12,157,105bp (see Table \ref{tbl:yeastchrs}).

Instead of copying all the folders and duplicate all of its files, we
will re-use from previous exercise by creating symbolic links
(something similar to a MS-Windows direct access) with the commands
below:

```{.sh}
#
# We are about to link previous exercise folders,
# assuming you have already set your current working directory
# on the exercise_04 folder...
#
for FD in seqs bowtie soapdenovo;
  do {
   ln -vs ../exercise_03/$FD ./;
  }; done

tail -6 $WDR/soapdenovo/SRR6130428/SRR6130428_k63_contig.log 
# There are 3156 contig(s) longer than 100, sum up 11779878 bp, with average length 3732.
# The longest length is 69946 bp, contig N50 is 14294 bp,contig N90 is 3279 bp.
# 4867 contig(s) longer than 64 output.
# Time spent on constructing contig: 0m.

# IMPORTANT: We have this data already from previous exercise,
#            we include here for completeness shake.
URL="https://downloads.yeastgenome.org/sequence/S288C_reference/genome_releases"
wget ${URL}/S288C_reference_genome_Current_Release.tgz \
     -O $WDR/seqs/S288C_reference_genome_Current_Release.tgz
     
pushd $WDR/seqs/
tar -zxvf S288C_reference_genome_Current_Release.tgz
popd

export GENOMEDIR=$WDR/seqs/S288C_reference_genome_R64-5-1_20240529;
export GENOMEFSA=S288C_reference_sequence_R64-5-1_20240529;
# you may consider copying those vars to your projectvars.sh file
               
zcat $GENOMEDIR/$GENOMEFSA.fsa.gz | \
  infoseq -only -length -noheading -sequence fasta::stdin  2> /dev/null | \
  gawk '{ s+=$1 } END{ print "# Yeast genome size (v.R64): "s }' 
# Yeast genome size (v.R64): 12157105
```

```{.sh}
#
###### NOTE ###### -------------------------------------------------------
#
# You can include a LaTeX table with the chromosome sizes and GC content,
# which can be produced with a command like this:
#
zcat $GENOMEDIR/$GENOMEFSA.fsa.gz | \
  infoseq -noheading -sequence fasta::stdin  2> /dev/null | \
    gawk 'BEGIN {
            printf "%15s & %10s & %12s & %10s \\\\\n",
                   "Chromosome", "GenBank ID", "Length (bp)", "GC content";
          }
          $0 != /^[ \t]*$/ {
	        L=$0;
	        sub(/^.*\[(chromosome|location)=/,"",L);
	        sub(/\].*$/,"",L);
	        sub(/_/,"\\_",$3);
            printf "%15s & %10s & %12d & %8.2f\\%% \\\\\n",
                   L, $3, $6, $7;
          }' > $WDR/docs/chromosomes_info.tex;
#
###### NOTE ###### -------------------------------------------------------
#
```

\newpage
```{=latex}
% if pandoc complains or it is not running LaTeX code as expected,
% this block can also be saved into a docs/fig_kmercountdist.tex file
% and then loaded into main doc using the LaTeX "\input" macro.
\begin{table}[!ht]
\begin{center}
\begin{tabular}{rrrr}
%%%
%%% Include your table here...
\csname @@input\endcsname docs/chromosomes_info
%%%
\end{tabular}
\parbox{0.75\linewidth}{%
 \caption[Reference \Scer\ chromosome summary.]{%
  \label{tbl:yeastchrs}\textbf{Reference \Scer\ chromosome summary.} This table will show information about length and GC content for the chromosomes of the \scer\ reference genome.
 }%caption
}%parbox
\end{center}
\end{table}
```

```{.sh}
# We need to compute some extra stats from the assemblies,
# we can download from github one of the assemblathon scripts for this purpose
# (to facilitate the task, it is already available on the bin folder)
#
## GITASM=https://github.com/KorfLab/Assemblathon/raw/master
## wget $GITASM/assemblathon_stats.pl  -O bin/assemblathon_stats_unfixed.pl
## wget http://korflab.ucdavis.edu/Unix_and_Perl/FAlite.pm -O bin/FAlite.pm
## chmod a+x bin/assemblathon_stats.pl
#
# however the main script fails due to a syntax error due to updates on Perl interpreter,
# here you can see the fixes introduced to the original script:
#
# 301c283
# < 	foreach my $size qw(1000 10000 100000 1000000 10000000){
# ---
# > 	foreach my $size (qw(1000 10000 100000 1000000 10000000)) {
# 428c409
# < 	foreach my $base qw (A C G T N){
# ---
# > 	foreach my $base (qw(A C G T N)) {

#
# IMPORTANT: run this first, otherwise you will get an error
#            about perl not finding FAlite.pm module
#
export PERL5LIB=$BIN;

#
zcat $GENOMEDIR/$GENOMEFSA.fsa.gz | \
  $BIN/assemblathon_stats.pl \
            -csv -genome_size 12160000 - \
          > $WDR/stats/assembly_stats_$GENOMEFSA.txt

$BIN/assemblathon_stats.pl \
            -csv -genome_size 12160000 \
            $WDR/soapdenovo/SRR6130428/SRR6130428_k63_graph.contig \
          > $WDR/stats/assembly_stats_SRR6130428_soapdenovo_k63_graph.contig.txt
          
# SRR6130428_k63_graph.contig.fa
```

We will continue using the smaller dataset from the previous exercise,
`SRR6130428`, for most of the examples on this exercise. However, if
you had already chosen another \sra\ reads set to generate your
assembly; you will need it to complete the table below with the
outcomes from the initial assembly steps (see the example from Table
\ref{tbl:assemblysummary}). If you ran assembly over a second set or
on a different reads sampling, you are also welcome to include such
information as another summary row on this table.

```{=latex}
% if pandoc complains or it is not running LaTeX code as expected,
% this block can also be saved into a docs/fig_kmercountdist.tex file
% and then loaded into main doc using the LaTeX "\input" macro.
\begin{table}[!ht]
\begin{center}
\begin{footnotesize}
\begin{tabular}{l|crrrrrrr}
%%%
%%% Your data here...
%%%
\bf Assembly          &  \bf Status & \bf\#Seqs & \bf Total Length & \bf Longest Seq & \bf Shortest Seq & \bf Avg Length & $N_{50}$ & $N_{50}$ \\
                      &             &           &           \it bp &          \it bp &           \it bp &         \it bp &  \#Seqs. &   \it bp \\  \hline \hline
Scer\_ref\_vR64       & Chromosomes &        17 &         12157105 &         1531933 &            85779 &         715124 &        6 &   924431 \\
SRR6130428            &     Contigs &      4867 &         11905871 &           69946 &               64 &            256 &    14236 &    14236 \\
\end{tabular}
\end{footnotesize}
\parbox{0.75\linewidth}{%
 \caption[Summary of assembly results on each set.]{%
  \label{tbl:assemblysummary}\textbf{Summary of assembly results on each set.} We have computed same parameters over \scer\ reference genome as for the assembled versions.
 }%caption
}%parbox
\end{center}
\end{table}
```

For some of the following analyses, we will focus on a couple of
reference genome chromosomes, `chrI` and `chrM` (mitochondrion
genome). Here we are going to filter them out:

```{.sh}
mkdir $WDR/seqs/chrs;

for SQ in chrI:NC_001133 chrM:NC_001224;
  do {
    SQN=${SQ%%:*}; # get the chr from SQ string
    SQI=${SQ##*:}; # get the refseq id from SQ string
    echo "# Filtering $SQN [$SQI] from whole genome fasta file..." 1>&2;
    echo 'ref|'$SQI'|' | \
      seqtk subseq $GENOMEDIR/$GENOMEFSA.fsa.gz - | \
        sed 's/^>.*$/>Scer_'$SQN'/;' \
                 > $WDR/seqs/chrs/$SQN.fa;
  }; done
```


# Exploring the assemblies

## Filter out contigs mapping to reference chromosomes

Later on we are going to run `dnadiff` from the `MUMmer`
package\footnote{S. Kurtz, A. Phillippy, A.L. Delcher, M. Smoot,
M. Shumway, C. Antonescu, and S.L. Salzberg.\newline\hspace*{1cm}
"Versatile and open software for comparing large genomes."
\textit{Genome Biology}, 5:R12, 2004.} to compare assembled contigs
against the reference or between them. In order to speed up the
example, we will first use `NCBI-BLAST`\footnote{C. Camacho,
G. Coulouris, V. Avagyan, N. Ma, J. Papadopoulos, K. Bealer, and
T.L. Madden.\newline\hspace*{1cm} "BLAST+: architecture and
applications." \textit{BMC Bioinformatics}, 10:421, 2008.} to project
all the assembled contigs into the chosen reference chromosomes, so we
will reduce all the downstream calculations. We need to create a
database for each assembly sequence sets and we will query by the
reference selected chromosomes.

```{.sh}
export SQSET=SRR6130428

mkdir -vp $WDR/blast/dbs

makeblastdb -in $WDR/soapdenovo/$SQSET/${SQSET}_k63_graph.contig.fa \
            -dbtype nucl \
	    -title "${SQSET}_SOAPdenovo_k63_contigs" \
	    -out $WDR/blast/dbs/${SQSET}_SOAPdenovo_k63_contigs \
	      2> $WDR/blast/dbs/${SQSET}_SOAPdenovo_k63_contigs.log 1>&2;
```

Now we can `BLAST` reference sequences against the newly created
database; we will use `megablast` option as we are comparing sequences
for the same species, and that `BLAST` program has the parameters
optimized for this kind of genomic searches.

```{.sh}
# we define here a custom BLAST tabular output format
BLASTOUTFORMAT='6 qseqid qlen sseqid slen qstart qend sstart send length';
BLASTOUTFORMAT=$BLASTOUTFORMAT' score evalue bitscore pident nident ppos positive';
BLASTOUTFORMAT=$BLASTOUTFORMAT' mismatch gapopen gaps qframe sframe';
export BLASTOUTFORMAT;

for SQ in chrI chrM;
  do {
    echo "# Running MEGABLAST: $SQSET x $SQ ..." 1>&2;
    mkdir -vp $WDR/blast/${SQ}-x-${SQSET};
    blastn -task megablast -num_threads 8                        \
           -db    $WDR/blast/dbs/${SQSET}_SOAPdenovo_k63_contigs \
           -outfmt "$BLASTOUTFORMAT"                             \
           -query $WDR/seqs/chrs/$SQ.fa                          \
           -out   $WDR/blast/${SQ}-x-${SQSET}/${SQ}-x-${SQSET}_megablast.out \
             2>   $WDR/blast/${SQ}-x-${SQSET}/${SQ}-x-${SQSET}_megablast.log;
  }; done
```

Once we have found the matching contigs, let's filter them out from
the whole assembly fasta file.

```{.sh}
for SQ in chrI chrM;
  do {
    echo "# Getting $SQSET contigs matching $SQ ..." 1>&2;
    OFBN="$WDR/blast/${SQ}-x-${SQSET}/${SQ}-x-${SQSET}_megablast";
    # this is a check to ensure we start with an empty contigs fasta file
    if [ -e "$OFBN.fa" ];
      then
        printf '' > "$OFBN.fa"; # rm can be also used here, but dangerous for novice
      fi;
    # get the contig IDs from third column and filter out sequences
    gawk '{ print $3 }' "$OFBN.out" | sed 's/^>//' | sort | uniq | \
      while read SQID;
        do {
           samtools faidx \
                    $WDR/soapdenovo/$SQSET/${SQSET}_k63_graph.contig.fa \
                    "$SQID";
        }; done >> "$OFBN.fa" \
                2> "$OFBN.fa.log";
  }; done
```

Now that we have the filtered out the required contigs, we can start
the following sequence comparison procedure based on `dnadiff`.

```{.sh}
for SQ in chrI chrM;
  do {
    #
    printf "# Running DNADIFF protocol: $SQSET x $SQ ..." 1>&2;
    IFBN="$WDR/blast/${SQ}-x-${SQSET}/${SQ}-x-${SQSET}_megablast";
    OFBD="$WDR/blast/${SQ}-x-${SQSET}/${SQ}-x-${SQSET}_megablast.dnadiff";
    #
    dnadiff -p $OFBD              \
            $WDR/seqs/chrs/$SQ.fa \
            $IFBN.fa              \
         2> $OFBD.log;
    printf " DNAdiff..." 1>&2;
    #
    mummerplot --large --layout --fat --postscript \
               -t "Alignment Plot: ${SQ}-x-${SQSET}" \
               -p $OFBD           \
                  $OFBD.1delta    \
               2> $OFBD.alnplot.log;
    # add this to convert PostScript image to PNG
    convert-im6 -verbose $OFBD.ps $OFBD.png;
    printf " ALNplot..." 1>&2;
    mummerplot --large --layout --fat --coverage --postscript \
               -t "Coverage Plot: ${SQ}-x-${SQSET}" \
               -p $OFBD.covg      \
                  $OFBD.1delta    \
               2> $OFBD.cvgplot.log;
    # add this to convert PostScript image to PNG
    convert-im6 -verbose $OFBD.covg.ps $OFBD.covg.png;
    printf " CVGplot..." 1>&2;
    #
    ( grep "TotalBases"   $OFBD.report;
      grep "AlignedBases" $OFBD.report;
      grep "AvgIdentity"  $OFBD.report
      ) > $OFBD.shortsummary;
    printf " DONE\n" 1>&2;
    #
  }; done
```

```{=latex}
\begin{figure}[!ht]
\begin{center}
 \begin{tabular}{c@{}c}
  \includegraphics[width=0.435\linewidth]{{blast/chrI-x-SRR6130428/chrI-x-SRR6130428_megablast.dnadiff}.png}     &
  \includegraphics[width=0.435\linewidth]{{blast/chrM-x-SRR6130428/chrM-x-SRR6130428_megablast.dnadiff}.png}     \\
  \includegraphics[width=0.435\linewidth]{{blast/chrI-x-SRR6130428/chrI-x-SRR6130428_megablast.dnadiff.covg}.png} &
  \includegraphics[width=0.435\linewidth]{{blast/chrM-x-SRR6130428/chrM-x-SRR6130428_megablast.dnadiff.covg}.png} \\
 \end{tabular}
 \parbox{0.8\linewidth}{%
  \caption[\texttt{dnadiff} comparison between two reference chromosomes and \texttt{SRR6130428} contigs]{%
   \label{fig:fastqcsetA}\textbf{\texttt{dnadiff} comparison between two reference chromosomes and \texttt{SRR6130428} contigs.} 
   Top panels show the alignment plots of contigs from \texttt{SRR6130428} assembly mapped over two reference \Scer\ chromosomes, 
   chrI and chrM on left and right panels respectively. Bottom panels show the alignment coverage for the same sequence relations. 
   It is evident from the comparison that contigs aligning to chrM have better contiguity despite overall coverages are quite 
   similar on both reference chromosomes.
  }%caption
 }%parbox
\end{center}
\end{figure}
``` 

Discuss later on which of the two chromosome assemblies do you think
had a better outcome and try to guess why. We may need to analyze
whether the contigs fully align to the reference chromosomes, __can you
calculate the average coverage of the aligned reads?__


## Assessment of genome completeness with `BUSCO`

`BUSCO`\footnote{M. Manni, M.R. Berkeley, M. Seppey, F.A. Simão, and
E.M. Zdobnov.\newline\hspace*{1cm} "\texttt{BUSCO} Update: Novel and
Streamlined Workflows along with Broader and Deeper Phylogenetic
Coverage for Scoring of Eukaryotic, Prokaryotic, and Viral Genomes."
\textit{Molecular Biology and Evolution}, 38(10)4647-–4654, 2021.}
estimates the completeness and redundancy of processed genomic data
based on universal single-copy orthologs. First of all, we need to
check if there is a clade-specific parameters set that fits with our
organism, \Scer\ belong to the _Saccharomycetes_ class, within the
_Ascomycota_ phylum in the Fungi kingdom (see
\href{https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=info&id=4932}{NCBI
Taxonomy browser species card}).

```{.sh}
busco --list-datasets
# 2022-10-24 19:06:49 INFO:	Downloading information on latest versions of BUSCO data...
# 2022-10-24 19:06:52 INFO:	Downloading file 'https://busco-data.ezlab.org/v5/data/information/lineages_list.2021-12-14.txt.tar.gz'
# 2022-10-24 19:06:53 INFO:	Decompressing file '/home/lopep/SANDBOX/BScCG2425/exercise_04/busco_downloads/information/lineages_list.2021-12-14.txt.tar.gz'
# 
# ################################################
# 
# Datasets available to be used with BUSCO v4 and v5:
# 
#  bacteria_odb10
#      - ...
#  archaea_odb10
#      - ...
#  eukaryota_odb10
#      - ...
#      - fungi_odb10
#          - ascomycota_odb10
#              - ...
#              - saccharomycetes_odb10
#              - ...
#          - ...
#      - ...
#  viruses (no root dataset)
#      - ...
# 
```

Luckily for us, there is one class specific set to evaluate the
completeness of our genome assembly `saccharomycetes_odb10`. However,
it is worth that you consider also to run the `BUSCO` commands below
using a more general parameters set, such as `fungi_odb10` or even
`eukaryota_odb10` (or both), and compare the corresponding
results. This can be useful to extrapolate the assessment to other
species genomes for which we will not have as much information as for
this yeast.

```{.sh}
mkdir -vp $WDR/busco/$SQSET

# fix PATH to point bin folders where you installed bbmap and busco programs
export PATH=$BIN:$BIN/bbmap:$BIN/busco/bin:$PATH

# this is to avoid a busco error complaining about paths in output file name
pushd $WDR/busco/$SQSET
# running busco in the output folder
busco -m genome \
      -i $WDR/soapdenovo/$SQSET/${SQSET}_k63_graph.contig.fa \
      -o ${SQSET}_k63_contigs \
      -l saccharomycetes_odb10 \
      -f;
# returning to the previous folder
popd;
##
## if you get an error "A run with the name xxx already exists"
## you should add command-line option "-f" to force overwriting those files
#
#    ---------------------------------------------------
#    |Results from dataset saccharomycetes_odb10       |
#    ---------------------------------------------------
#    |  C:98.8%[S:96.7%,D:2.1%],F:0.5%,M:0.7%,n:2137   |
#    |2111  Complete BUSCOs (C)                        |
#    |2067  Complete and single-copy BUSCOs (S)        |
#    |44    Complete and duplicated BUSCOs (D)         |
#    |11    Fragmented BUSCOs (F)                      |
#    |15    Missing BUSCOs (M)                         |
#    |2137  Total BUSCO groups searched                |
#    ---------------------------------------------------
# 2022-10-26 19:55:26 INFO: BUSCO analysis done. Total running time: 440 seconds
```

You can take a step further to plot the resulting completeness values,
using the `generate_plot.py` script that is provided under the scripts
folder of your `busco` installation:

```{.sh}
# If you installed BUSCO from sources, this script should be at:
#    $BIN/busco/scripts/generate_plot.py
# however if you got from the conda/mamba environment you must
# search on the conda env bin folder; check your installation
# for the equivalent file, i.e.: 
#    $HOME/.conda/envs/BScBI-CG2425_exercises/bin/generate_plot.py
python3 /home/maria/busco/scripts/generate_plot.py \
    -wd $WDR/busco/$SQSET/${SQSET}_k63_contigs
```

**Embed here the corresponding figure**; if you have ran more than one
assembly, you can combine the `busco` results to compare completeness
among them.

```{=latex}
% you can also use an \input macro to load the figure from a tex file.
%
%%% YOUR FIGURE HERE %%%

\begin{figure}[!ht]
\begin{center}
\includegraphics[width=0.8\linewidth]{{busco/SRR6130428/SRR6130428_k63_contigs/busco_figure}.png}
\parbox{0.75\linewidth}{%
 \caption[BUSCO Analysis]{%
  \label{fig:fastqcsetA}\textbf{BUSCO Analysis.}
  Plot for the resulting completeness values of the BUSCO ANALYSIS. Shows a high degree of completeness
  of complete contigs (98.9%), fragmented contigs (0.5%), and percentage of missing coontigs (0.7%). 
 }%caption
}%parbox
\end{center}
\end{figure}
```
There is a high degree
of completeness, since the percentage of complete contigs is 98.8%, the percentage
of fragmented contigs is 0.5% and the percentage of missing contigs is 0.7%.

# Masking assembly repetitive sequences

Although installing `RepeatMasker`\footnote{A.F.A. Smit, R. Hubley,
and P. Green, "RepeatMasker at \url{http://repeatmasker.org}"
(\textit{unpublished}).} could be a complex task, specially if we want
to use the `RepBase` database of curated repeats, we can just install
one of the search engines (`hmmer` in our case for this exercise). We
are going to run `RepeatMasker` on the set of contigs that map to the
reference chromosome chrM for testing purposes (looking that
everything is working); we are going to use some extra parameters to
speed up the searches (\texttt{-qq}, much faster but 10\% less
sensitive), to use lowercase nucleotides to mask repeats on the output
masked sequences (\texttt{-small}), and to get extra output files
(\texttt{-gff}, for the repeats location coords in GFF).

```{.sh}
SQ=chrM;
SQSET=SRR6130428;
CONDABIN=$HOME/miniconda3/envs/BScBI-CG2425_exercises/bin;
mkdir -vp $WDR/repeatmasker;

# Run first on the chrM contigs set for a quick test
#
mkdir -vp $WDR/repeatmasker/${SQ}-x-${SQSET}_megablast
$CONDABIN/RepeatMasker \
           $WDR/blast/${SQ}-x-${SQSET}/${SQ}-x-${SQSET}_megablast.fa \
           -qq -small -gff -species 'Saccharomyces' -e hmmer \
           -dir $WDR/repeatmasker/${SQ}-x-${SQSET}_megablast    \
             2> $WDR/repeatmasker/${SQ}-x-${SQSET}_megablast.rptmskr.log 1>&2

# You can run now over the whole genome assemblies # CAUTION! it may take too long
#
mkdir -vp $WDR/repeatmasker/${SQSET}_k63_graph.contig
$CONDABIN/RepeatMasker \
           $WDR/soapdenovo/$SQSET/${SQSET}_k63_graph.contig.fa \
           -qq -small -gff -species 'Saccharomyces' -e hmmer \
           -dir $WDR/repeatmasker/${SQSET}_k63_graph.contig \
             2> $WDR/repeatmasker/${SQSET}_k63_graph.contig.rptmskr.log 1>&2

# You can also mask the second assembly...

SQ=chrI;
SQSET=SRR6130428;
CONDABIN=$HOME/miniconda3/envs/BScBI-CG2425_exercises/bin;
mkdir -vp $WDR/repeatmasker;

# Run first on the chrM contigs set for a quick test
#
mkdir -vp $WDR/repeatmasker/${SQ}-x-${SQSET}_megablast
$CONDABIN/RepeatMasker \
           $WDR/blast/${SQ}-x-${SQSET}/${SQ}-x-${SQSET}_megablast.fa \
           -qq -small -gff -species 'Saccharomyces' -e hmmer \
           -dir $WDR/repeatmasker/${SQ}-x-${SQSET}_megablast    \
             2> $WDR/repeatmasker/${SQ}-x-${SQSET}_megablast.rptmskr.log 1>&2

# You can run now over the whole genome assemblies # CAUTION! it may take too long
#
mkdir -vp $WDR/repeatmasker/${SQSET}_k63_graph.contig
$CONDABIN/RepeatMasker \
           $WDR/soapdenovo/$SQSET/${SQSET}_k63_graph.contig.fa \
           -qq -small -gff -species 'Saccharomyces' -e hmmer \
           -dir $WDR/repeatmasker/${SQSET}_k63_graph.contig \
             2> $WDR/repeatmasker/${SQSET}_k63_graph.contig.rptmskr.log 1>&2
```

**Include a summary table here** listing the repeats families found with
the total counts. Describe the number of masked nucleotides too.

```{=latex}
% you can also use an \input macro to load the table from a tex file.
%
%%% YOUR TABLE/RESULTS HERE %%%

\begin{table}
 \begin{center}
  \begin{small} % to get a small font and fit table within the margins...
   \begin{tabular}{lcrr}
   & & \textbf{Length} & \textbf{Percentage of Sequence} \\ \hline\hline
   \textbf{File Name:} & \multicolumn{3}{l}{\texttt{chrM-x-SRR6130428\_megablast.fa}} \\
   \textbf{Sequences:} & 45 & & \\
   \textbf{Total Length:} & 119699 bp & & \\
   \textbf{GC Level:} & 22.83\% & & \\
   \textbf{Bases Masked:} & 9879 bp & (8.25\%) & \\ \hline
   \hline
   \textbf{Category} & \textbf{Number of Elements} & \textbf{Length Occupied} & \textbf{Percentage of Sequence} \\ \hline\hline
   Retroelements & 0 & 0 bp & 0.00\% \\
   \quad SINEs & 0 & 0 bp & 0.00\% \\
   \quad Penelope & 0 & 0 bp & 0.00\% \\
   \quad LINEs & 0 & 0 bp & 0.00\% \\
   \quad\quad CRE/SLACS & 0 & 0 bp & 0.00\% \\
   \quad\quad L2/CR1/Rex & 0 & 0 bp & 0.00\% \\
   \quad\quad R1/LOA/Jockey & 0 & 0 bp & 0.00\% \\
   \quad\quad R2/R4/NeSL & 0 & 0 bp & 0.00\% \\
   \quad\quad RTE/Bov-B & 0 & 0 bp & 0.00\% \\
   \quad\quad L1/CIN4 & 0 & 0 bp & 0.00\% \\
   \quad LTR elements & 0 & 0 bp & 0.00\% \\
   \quad\quad BEL/Pao & 0 & 0 bp & 0.00\% \\
   \quad\quad Ty1/Copia & 0 & 0 bp & 0.00\% \\
   \quad\quad Gypsy/DIRS1 & 0 & 0 bp & 0.00\% \\
   \quad\quad\quad Retroviral & 0 & 0 bp & 0.00\% \\ \hline
   DNA Transposons & 0 & 0 bp & 0.00\% \\
   \quad hobo-Activator & 0 & 0 bp & 0.00\% \\
   \quad Tc1-IS630-Pogo & 0 & 0 bp & 0.00\% \\
   \quad En-Spm & 0 & 0 bp & 0.00\% \\
   \quad MULE-MuDR & 0 & 0 bp & 0.00\% \\
   \quad PiggyBac & 0 & 0 bp & 0.00\% \\
   \quad Tourist/Harbinger & 0 & 0 bp & 0.00\% \\
   \quad Other (Mirage, P-element, Transib) & 0 & 0 bp & 0.00\% \\ \hline
   Rolling-circles & 0 & 0 bp & 0.00\% \\ \hline
   Unclassified & 0 & 0 bp & 0.00\% \\ \hline
   Total interspersed repeats & & 0 bp & 0.00\% \\ \hline
   Small RNA & 0 & 0 bp & 0.00\% \\ \hline
   Satellites & 0 & 0 bp & 0.00\% \\
   Simple repeats & 148 & 8837 bp & 7.38\% \\
   Low complexity & 11 & 1042 bp & 0.87\% \\ \hline
  \end{tabular}
  \end{small}
 \end{center} 
 \caption[Summary of sequence composition for \texttt{chrM-x-SRR6130428\_megablast.fa}]{%
 \label{tbl:sequence_composition}\textbf{Summary of sequence composition for file \texttt{chrM-x-SRR6130428\_megablast.fa}.} 
 Overview of various types of repeats, transposons, and base composition across the sequence, with details on their lengths and percentages.
}
\end{table}
```

```{=latex}
\begin{table}
 \begin{center}
  \begin{small} % to get a small font and fit table within the margins...
   \begin{tabular}{lcrr}
   & & \textbf{Length} & \textbf{Percentage of Sequence} \\ \hline\hline
   \textbf{File Name:} & \multicolumn{3}{l}{\texttt{chrI-x-SRR6130428\_megablast.fa}} \\
   \textbf{Sequences:} & 500 & & \\
   \textbf{Total Length:} & 797359 bp & & \\
   \textbf{GC Level:} & 38.19\% & & \\
   \textbf{Bases Masked:} & 10534 bp & (1.32\%) & \\ \hline
   \hline
   \textbf{Category} & \textbf{Number of Elements} & \textbf{Length Occupied} & \textbf{Percentage of Sequence} \\ \hline\hline
   Retroelements & 0 & 0 bp & 0.00\% \\
   \quad SINEs & 0 & 0 bp & 0.00\% \\
   \quad Penelope & 0 & 0 bp & 0.00\% \\
   \quad LINEs & 0 & 0 bp & 0.00\% \\
   \quad\quad CRE/SLACS & 0 & 0 bp & 0.00\% \\
   \quad\quad L2/CR1/Rex & 0 & 0 bp & 0.00\% \\
   \quad\quad R1/LOA/Jockey & 0 & 0 bp & 0.00\% \\
   \quad\quad R2/R4/NeSL & 0 & 0 bp & 0.00\% \\
   \quad\quad RTE/Bov-B & 0 & 0 bp & 0.00\% \\
   \quad\quad L1/CIN4 & 0 & 0 bp & 0.00\% \\
   \quad LTR elements & 0 & 0 bp & 0.00\% \\
   \quad\quad BEL/Pao & 0 & 0 bp & 0.00\% \\
   \quad\quad Ty1/Copia & 0 & 0 bp & 0.00\% \\
   \quad\quad Gypsy/DIRS1 & 0 & 0 bp & 0.00\% \\
   \quad\quad\quad Retroviral & 0 & 0 bp & 0.00\% \\ \hline
   DNA Transposons & 0 & 0 bp & 0.00\% \\
   \quad hobo-Activator & 0 & 0 bp & 0.00\% \\
   \quad Tc1-IS630-Pogo & 0 & 0 bp & 0.00\% \\
   \quad En-Spm & 0 & 0 bp & 0.00\% \\
   \quad MULE-MuDR & 0 & 0 bp & 0.00\% \\
   \quad PiggyBac & 0 & 0 bp & 0.00\% \\
   \quad Tourist/Harbinger & 0 & 0 bp & 0.00\% \\
   \quad Other (Mirage, P-element, Transib) & 0 & 0 bp & 0.00\% \\ \hline
   Rolling-circles & 0 & 0 bp & 0.00\% \\ \hline
   Unclassified & 0 & 0 bp & 0.00\% \\ \hline
   Total interspersed repeats & & 0 bp & 0.00\% \\ \hline
   Small RNA & 0 & 0 bp & 0.00\% \\ \hline
   Satellites & 0 & 0 bp & 0.00\% \\
   Simple repeats & 206 & 9022 bp & 1.13\% \\
   Low complexity & 32 & 1512 bp & 0.19\% \\ \hline
  \end{tabular}
  \end{small}
 \end{center} 
 \caption[Summary of sequence composition for \texttt{chrI-x-SRR6130428\_megablast.fa}]{%
 \label{tbl:sequence_composition_chrI}\textbf{Summary of sequence composition for file \texttt{chrI-x-SRR6130428\_megablast.fa}.} 
 Overview of various types of repeats, transposons, and base composition across the sequence, with details on their lengths and percentages.
}
\end{table}
```

```{=latex}
\begin{table}
 \begin{center}
  \begin{small} % to get a small font and fit table within the margins...
   \begin{tabular}{lcrr}
   & & \textbf{Length} & \textbf{Percentage of Sequence} \\ \hline\hline
   \textbf{File Name:} & \multicolumn{3}{l}{\texttt{SRR6130428\_k63\_graph.contig.fa}} \\
   \textbf{Sequences:} & 4867 & & \\
   \textbf{Total Length:} & 11905871 bp & & \\
   \textbf{GC Level:} & 38.12\% & & \\
   \textbf{Bases Masked:} & 159289 bp & (1.34\%) & \\ \hline
   \hline
   \textbf{Category} & \textbf{Number of Elements} & \textbf{Length Occupied} & \textbf{Percentage of Sequence} \\ \hline\hline
   Retroelements & 0 & 0 bp & 0.00\% \\
   \quad SINEs & 0 & 0 bp & 0.00\% \\
   \quad Penelope & 0 & 0 bp & 0.00\% \\
   \quad LINEs & 0 & 0 bp & 0.00\% \\
   \quad\quad CRE/SLACS & 0 & 0 bp & 0.00\% \\
   \quad\quad L2/CR1/Rex & 0 & 0 bp & 0.00\% \\
   \quad\quad R1/LOA/Jockey & 0 & 0 bp & 0.00\% \\
   \quad\quad R2/R4/NeSL & 0 & 0 bp & 0.00\% \\
   \quad\quad RTE/Bov-B & 0 & 0 bp & 0.00\% \\
   \quad\quad L1/CIN4 & 0 & 0 bp & 0.00\% \\
   \quad LTR elements & 0 & 0 bp & 0.00\% \\
   \quad\quad BEL/Pao & 0 & 0 bp & 0.00\% \\
   \quad\quad Ty1/Copia & 0 & 0 bp & 0.00\% \\
   \quad\quad Gypsy/DIRS1 & 0 & 0 bp & 0.00\% \\
   \quad\quad\quad Retroviral & 0 & 0 bp & 0.00\% \\ \hline
   DNA Transposons & 0 & 0 bp & 0.00\% \\
   \quad hobo-Activator & 0 & 0 bp & 0.00\% \\
   \quad Tc1-IS630-Pogo & 0 & 0 bp & 0.00\% \\
   \quad En-Spm & 0 & 0 bp & 0.00\% \\
   \quad MULE-MuDR & 0 & 0 bp & 0.00\% \\
   \quad PiggyBac & 0 & 0 bp & 0.00\% \\
   \quad Tourist/Harbinger & 0 & 0 bp & 0.00\% \\
   \quad Other (Mirage, P-element, Transib) & 0 & 0 bp & 0.00\% \\ \hline
   Rolling-circles & 0 & 0 bp & 0.00\% \\ \hline
   Unclassified & 0 & 0 bp & 0.00\% \\ \hline
   Total interspersed repeats & & 0 bp & 0.00\% \\ \hline
   Small RNA & 0 & 0 bp & 0.00\% \\ \hline
   Satellites & 0 & 0 bp & 0.00\% \\
   Simple repeats & 3248 & 135605 bp & 1.14\% \\
   Low complexity & 505 & 23684 bp & 0.20\% \\ \hline
  \end{tabular}
  \end{small}
 \end{center} 
 \caption[Summary of sequence composition for \texttt{SRR6130428\_k63\_graph.contig.fa}]{%
 \label{tbl:sequence_composition_k63}\textbf{Summary of sequence composition for file \texttt{SRR6130428\_k63\_graph.contig.fa}.} 
 Overview of various types of repeats, transposons, and base composition across the sequence, with details on their lengths and percentages.
}
\end{table}
```

\newpage
# Discussion

__IMPORTANT__ Discuss your results here (around 300 words). And
remember to include in the Appendices section (see page
\pageref{sec:appendices}), any extra script you wrote from this
exercise `bin` folder using the `loadfile` macro. We can take
advantage of the \LaTeX\ referencing capabilities, as described in the
first exercise template.

To start our analysis, we first explored the assemblies and compare the assembled contigs against the reference or between
them. To do so, we projected the assembled contigs into the chosen reference chromosomes to reduce the downstream calculations. 
Then, we performed a BLAST between the reference sequence and a newly created database, found the matching contigs, filtered
them and compared the sequences based on the procedure of dnadiff. In our resulting plots, we can see that the contigs aligning 
to chrM have a better contiguity than chrI despite both of the coverages being quite similar on both reference chromosomes. 

Furthermore, we continued by performing the BUSCO analysis to estimate the completeness and redundancy of processed genomic
data based on universal single-copy orthologs. So, we look at our plot on the resulting completeness, and we observe a 98.8% of 
complete contigs, a percentage of 0.5% of fragmented contigs and a 0.7% of missing contigs. 

Finally, we run RepeatMasker on a set of contigs that map to the reference chromosomes (chrI and chrM). We obtain 
three tables and these show the comparison of DNA sequence features across files `chrM`, `chrI`, and `SRR6130428_k63`. `chrM` has a lower GC 
content (22.83%) and a higher proportion of simple repeats (7.38%) compared to `chrI` (38.19% GC) and `SRR6130428_k63` 
(38.12% GC), suggesting `chrM` may represent mitochondrial DNA, which typically has lower GC and higher repeat content. 
In contrast, `chrI` and `SRR6130428_k63` are likely nuclear sequences, given their higher GC content and larger sequence counts.
Moreover, none of the files contain significant transposable elements like LINEs, SINEs, or LTRs, which may indicate regions with 
fewer mobile elements or pre-processed data with repetitive elements excluded. Despite varying sequence lengths and complexities,
`chrI` and `SRR6130428_k63` share similar levels of simple repeats and low complexity regions, consistent with typical 
nuclear genome structure. These patterns likely reflect the unique structural characteristics of mitochondrial vs. nuclear genomes.
To conclude, based on these results we can say that our genome assembly has high quality, with few errors and well documented
repetitive content. 

\clearpage

# Appendices
\label{sec:appendices}

## Software

We have used the following versions:

```{.sh}
uname -a
# Linux aleph 5.15.0-48-generic #54-Ubuntu SMP 
# Fri Aug 26 13:26:29 UTC 2022 x86_64 x86_64 x86_64 GNU/Linux

R --version
# R version 4.3.1 (2023-06-16) -- "Beagle Scouts"
# Copyright (C) 2023 The R Foundation for Statistical Computing
# Platform: x86_64-conda-linux-gnu (64-bit)

infoseq -version
# EMBOSS:6.6.0.0

wget --version
# GNU Wget 1.21.2 built on linux-gnu.

pandoc --version
# pandoc 3.1.3
# Features: +server +lua
# Scripting engine: Lua 5.4

mamba --version
# mamba 1.4.2
# conda 23.3.1

bunzip2 --version
# bzip2, a block-sorting file compressor.  Version 1.0.8, 13-Jul-2019.

mummerplot -v
# mummerplot 3.5

busco -v
# BUSCO 5.5.0

RepeatMasker -v
# RepeatMasker version 4.1.5
```


## Supplementary files
\label{sec:supplfiles}


### `conda` environment dependencies for the exercise

\loadfile{environment.yml}{environment.yml}{prg:environmentYML}


### Project specific scripts

```{=latex}
% now it is commented on the LaTeX document
% \loadfile{an\_script\_example.pl}{bin/an_script_example.pl}{prg:scriptexamplePERL}
```

### Shell global vars and settings for this project

\loadfile{projectvars.sh}{projectvars.sh}{prg:projectvarsBASH}


## About this document

This document was be compiled into a PDF using `pandoc` (see
`projectvars.sh` from previous subsection) and some `LaTeX` packages
installed in this linux system. `synaptic`, `apt-get` or `aptitude`
can be used to retrieve and install those tools from linux
repositories. As the `raw_tex` extension has been provided to the
`markdown_github` and `tex_math_dollars` formats, now this document
supports inline \LaTeX\ and inline formulas!

You can get further information from the following links about the
[Mark Down syntax](http://daringfireball.net/projects/markdown/syntax#link),
as well as from the manual pages (just type `man pandoc` and/or `man
pandoc_markdown`).
