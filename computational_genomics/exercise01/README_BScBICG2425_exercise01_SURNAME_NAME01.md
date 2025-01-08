---
### BEGIN-Of-YAML-Block ###
#
## ######################################################################################
##
##   README_BScBICG2425_exercise01_SURNAME_NAME.md
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
thyexercise: 01
#
# title-meta is title string without any LaTeX format to be used as pdftitle, part of emails subject...
title-meta: BScBI-CG Exercise 01
#
# title is the big title for the cover page, fully LaTeX formated to fit into a shortstack command...
title: |
  \textsc{BScBI-CG}
  \vskip 1.25ex \textsc{Practicals}
  \vskip 1.25ex \textsc{Report}
#subtitle:
#
# runtitle is the running header or footer, used i.e. by fancyheadings...
runtitle: |
  BScBI-CG Exercise 01 Report
#
# author-meta sets the pdfauthor variable...
author-meta: !!str 'Maria Cobo @ BScBI Computational Genomics'
#
# authors to appear on the title page...
author:
- name: Maria Cobo
  myemail: maria.cobo
  mydomain: alum.esci.upf.edu
# you can follow the example below
# we need to define email as two fields (user/domain) to avoid parsing email by YAML
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
-  \def\Dmel{\textit{Drosophila melanogaster}}
-  \def\dmel{\textit{D.~melanogaster}}
-  \def\Amel{\textit{Apis mellifera}}
-  \def\amel{\textit{A.~mellifera}}
#
### End-Of-YAML-Block ###
---

```{=comment}
We do not need the comment LaTeX environment to hide a block,
with pandoc >2.x one can define "comment" blocks, which are even
escaped and not passed to LaTeX and pandoc (achieving the same hidden effect).
\begin{comment}
The \texttt{comment} \LaTeX\ environment allow us to include any piece 
of text, code, etc, but it wil not be included into the final PDF report 
when compiling the \texttt{MarkDown} file. You can open/close one 
of such environments at any time if you need them.
\end{comment}
```

# Introduction

We are going to reuse an existing `MarkDown` report template, which
for this practical exercise describes a basic sequence analysis
protocol to explore a fly transcriptome dataset. Some parameters will
be estimated from the set of sequences provided. Following the same
protocol, you can provide a similar analysis for another species
trancriptome, and then compare the results on the discussion section
(see page \pageref{sec:discussion}).


## Objectives

* We will practice again how to report the commands and the results
  required to solve the practical exercises, using a text file, this
  report of course, and the `MarkDown` syntax. Some enhancements using
  \LaTeX\ will be introduced at every practical session.
* The main goal though is to reproduce a simple sequence analysis
  protocol, on the provided species transcriptome, and then extend to
  a second species transcriptome. This will make possible to discuss
  results separately but also to compare between each set.
* We will install and use basic programs from the Bioinformatics
  software suite `EMBOSS`, to calculate for instace GC content from
  the command-line interface.



## Prerequisites

### Installing required software

As for the previous practical, we must ensure that at least `pandoc`
and `pdflatex` commands are running smoothly over our report
files. You probably may need to install the following packages before
working on anything else, unless you got them running from the
previous exercise of course. We will need to install the `EMBOSS`
suite, but those instructions are described on the data retrieval
section (see section \pageref{sec:emboss} on page
\pageref{sec:emboss}). Just consider the different package manager
tools that are available on each system distribution:

```{.sh}
#################################
# Examples of package managers installing "emboss"
# on different operative systems, use only the ones suited for your operative system...

## For Linux distributions:

# on a debian/ubuntu/mint linux system (DEBs)
apt-cache search emboss     # to check if there is such a package
sudo apt-get install emboss # to install such a package

# on a redhat/fedora/centos linux system (RPMs)
yum search emboss           # to check if there is such a package
su -c 'yum install emboss'

# on a SUSE/openSuse  linux system
zypper search "emboss"
sudo zypper install emboss

## For MacOS distributions:

# on a Mac system using homebrew packages (**recommended option on a Mac**,
# see tutorial on the course introduction section materials at virtual campus)
brew search emboss
# check the above command output, i.e. "brewsci/bio/emboss", to use on install:
sudo brew install brewsci/bio/emboss

# on a Mac system using anaconda packages (https://conda.io/docs/index.html)
conda search emboss
# check the above command output to use on install:
sudo conda install -c bioconda emboss

# on a Mac system using mac ports (https://guide.macports.org/)
port search emboss
# check the above command output to use on install:
sudo port install emboss

## IMPORTANT ## Do not mess your Mac system using all
#               of the previous three install options, use the one
#               already available on your system or install "homebrew".

## Other:

# you can also install the package if available for the CygWin environment
# running on a Windows box (hhtp://www.cygwin.com/)

# add your packaging system here if you have not used any of the above commands...

# ...

#
#################################
```

#### __Linux__ users:

You must submit an exercise report for each practical as two single
files, a `MarkDown` text file and a `PDF` compiled from that
`MarkDown` file. In order to run such procedure, we must ensure that
we have the software tools and the corresponding dependencies. This
has to be done for the first exercise only, as you will have those
tools already available for future compilations of the other
exercises. The command below will install those tools:

```{.sh}
sudo apt-get install pandoc \
                     texlive-latex-recommended \
                     texlive-latex-extra \
                     texlive-fonts-recommended

# getting the latest version from repository
#     https://github.com/jgm/pandoc/releases/latest
# for a Debian-based linux distribution you can run those two commands:
wget https://github.com/jgm/pandoc/releases/download/3.1.8/pandoc-3.4-1-amd64.deb
sudo dpkg -i pandoc-3.4-1-amd64.deb
```

When using `pandoc` version greater than `2.x` we will be able to
apply further macros and tags on our report file (like embed `LaTeX`
blocks). In case you want to play with \LaTeX\, I will recommend you
to install the complete set with this command:

```{.sh}
sudo apt-get install texlive-full \
                     texlive-fonts-recommended \
                     texlive-fonts-extra
```

You can also install optional packages, such a text editor with
programming facilities and extensions, like `emacs` or `geany` (you
can also use `sublime`, `atom`, `gedit`, ...):

```{.sh}
sudo apt-get install emacs geany vim vim-gtk
```

Further instructions will be given on the templates in case a
practical requires that you install further software...

#### __MacOS__ users:

MacOS users can try to install ports for the software tools from
`homebrew` or `conda` repositories. Here we have a brief summary of
such commands. `Homebrew` is a package manager for __MacOS__; it also
has a port for __Linux__, known as `Linuxbrew`. To install this
package manager, just copy the following command and paste it to a
Terminal window:

```{.sh}
# setting up brew tool and its repositories
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install.sh)"
```

A list of all the available packages (known as *formulae*) is found at
[formulae.brew.sh](https://formulae.brew.sh). Examples of `brew`
command to install the `EMBOSS` or the `pandoc` are shown here:

```{.sh}
# Getting EMBOSS installed with brew
brew search emboss
# homebrew/science/emboss  #<-- check the output of the previous command to use in your system
brew install homebrew/science/emboss

# Getting pandoc installed with brew
brew install pandoc
# this is optional
brew install pandoc-citeproc
# you need this to typeset PDFs with LaTeX
brew install librsvg python homebrew/cask/basictex

# You can also download the pandoc MacOS instaler from the following URL:
# https://github.com/jgm/pandoc/releases/download/3.4/pandoc-3.4-x86_64-macOS.pkg
```

#### Using `conda/mamba` environments:

As we saw in the previous exercise, another way to install the
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
\pageref{prg:environmentYML}),

### Initializing the main report files

In a `bash` command-line terminal, create a folder for all the
practicals on this subject and change working dir to that folder:

```{.sh}
mkdir practicals
cd practicals
```

Then, we need to download the exercise tarball (the `*.tgz` file) from
the \CGVC\ into that folder, unpack such file, modify the files
accordingly to the user within the exercise folder, and set it as the
current working directory for the rest of the practical session...

```{.sh}
# Uncompress and unpack the exercise files from tarball
#
# NOTE: If you are reading this, you probably have already done this step.
#
tar -zxvf BScBI_CG2425_exercise_01.tgz

# Move into the new extracted folder
cd exercise_01

# Rename the MarkDown README_*.md file by replacing NAME and SURNAME strings
# with your "NAME" and "SURNAME" with the following command
mv -v README_BScBICG2425_exercise01_SURNAME_NAME.md \
      README_BScBICG2425_exercise01_YourSurname_YourName.md

# Open exercise files using your text editor of choice
# (for instance vim, emacs, gedit, sublime, atom, ...);
emacs  projectvars.sh \
       README_BScBICG2425_exercise01_YourSurname_YourName.md

# Fix "NAME" and "SURNAME" placeholders on those files
# and save the changes before continuing.

# Load the bash definitions from projectvars.sh
source projectvars.sh
# for instance the variable WDR was set to the absolute path
# to current exercise working directory
echo $WDR

# Now you are ready to start the practical by looking at
# the MarkDown file for further instructions to run the corresponding code blocks.

# Each time you include your answers/code/results on the README file, you can compile it into PDF.
# So that, let's tests if we can compile the modified MarkDown document.
# You probably must install some dependencies yet...
runpandoc
```

Again, remember to submit both files, the MD and the PDF, to the
\CGVC\, once errors/warnings have been fixed, all the requested task
have been completed, and you have discussed your results on the
corresponding sections.

If you have succeeded on the software installation step, then you can
start with the analyses provided on the next section... May the shell
be with you...

\newpage


# Genomic Analyses on Command-line

## Datasets

We have to analyze sequence length distribution and GC content for the
current mRNA sequences annotated on \Dmel\ genome (BDGP Release `6 +
ISO1 MT/dm6`, Aug. 2014). We have connected to the 
[UCSC genome browser download web site](http://hgdownload.soe.ucsc.edu/downloads.html),
and followed the
["Genome sequence files and select annotations (2bit, GTF, GC-content,
etc)" link](https://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/). 
We can see that there is a file named `mrna.fa.gz`, we can just copy
the corresponding link to the following command:

```{.sh}
export DT=$WDR/data
mkdir -v $DT
# You can also add the previous var definition to your 'projectvars.sh' file
# so it will be saved and can be easily reused when sourcing the file again.

wget https://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/mrna.fa.gz \
     -O $DT/dmel_mrna.fa.gz
```

We obtained a compressed file containing all the mRNA sequences
annotated for this version of the fly genome in `fasta` format. We can
check it's content:

```{.sh}
zcat $DT/dmel_mrna.fa.gz | head -10
# >DQ327735 1
# ttgttgcgcacacgcaccagaagagaggaggatcgaccaggcagctgttc
# tctggctctctggaaagtggtcaaaggagaaggaggaggtcgtaagagcg
# gtagaatcgacaatataatcggagtcatatcgggcatcaacgtcggccac
# atcaacatcaacaacggcagcagacgtcgctaattgcaaccaacacggga
# gctgcagcctggaccctacatatccaatgttcagaatttaaatgcaccac
# agcagcaacatcgaggtcctcagcagcagcatcagacggccggttgcctc
# aaagttccttagcctgcgtggcacgtggcatccatcagcgcctcaggctg
# aataaaaggaaatcccagccataaaaagtaaatcacttcggcgccaacca
# gttccagtcatcatggctcactccaaggtgatccttaagatcctgttttc
```

Let's see how many mRNA sequences do we have. We can uncompress the
file and keep going with the flat text file, which can take a lot of
disk resources, or we can keep using the compressed file as in the
previous command.

```{.sh}
zcat $DT/dmel_mrna.fa.gz | egrep -c '^>' 
# 75218
```

In case that a program can deal with compressed files, we can take
adantage of that feature. There is an `egrep` version that can read
such files, if you were guessing is `zegrep` of course (as it happens
with `cat` and `zcat`).

```{.sh}
zegrep -c '^>' $DT/dmel_mrna.fa.gz
# 75218
```


## Simple analysis of mRNA sequences
\label{sec:emboss}

To calculate the nucleotide length and the GC content, we can use one
of the programs that is provided within the [EMBOSS suite](http://emboss.sourceforge.net/):
`infoseq`. You can get information about this tool 
[from this link](http://emboss.sourceforge.net/apps/release/6.6/emboss/apps/infoseq.html).

```{.sh}
# You get info about this tool command-line options
infoseq -help
# If 'emboss-doc' has been installed youcan also use for further details
man infoseq

# Let's run it on our dataset
zcat $DT/dmel_mrna.fa.gz | \
  infoseq -sequence fasta::stdin \
          -noheading -only -name -length -pgc | \
    head -5
# Display basic information about sequences
# DQ327735       1883   53.48
# DQ327736       272    55.15
# DQ327737       338    60.36
# DQ327738       1294   50.85
# DQ327739       758    50.92
```

In the above command we are just looking to the expected output, the
following is doing the job and saving the output into
`dmel_mrna.lengc.tbl` file.

```{.sh}
zcat $DT/dmel_mrna.fa.gz | \
  infoseq -sequence fasta::stdin \
          -outfile $DT/dmel_mrna.lengc.tbl \
          -noheading -only -name -length -pgc
```

However, it is advisable to work as much as possible with compressed
files, so the following commands could be also useful.

```{.sh}
zcat $DT/dmel_mrna.fa.gz | \
  infoseq -sequence fasta::stdin \
          -noheading -only -name -length -pgc | \
    gzip -vc9 - > $DT/dmel_mrna.lengc.tbl.gz
```


## Visualizing the analysis

By running `R` command, we enter in the `R` shell interpreter, which
understands `R` commands of course.

```{.sh}
R
# 
# R version 4.1.2 (2021-11-01) -- "Bird Hippie"
# Copyright (C) 2021 The R Foundation for Statistical Computing
# Platform: x86_64-pc-linux-gnu (64-bit)
# 
# R is free software and comes with ABSOLUTELY NO WARRANTY.
# You are welcome to redistribute it under certain conditions.
# Type 'license()' or 'licence()' for distribution details.
# 
# R is a collaborative project with many contributors.
# Type 'contributors()' for more information and
# 'citation()' on how to cite R or R packages in publications.
# 
# Type 'demo()' for some demos, 'help()' for on-line help, or
# 'help.start()' for an HTML browser interface to help.
# Type 'q()' to quit R.
#
```

Now, we must load the tabular data into a variable.


```{.r}
# if we have an uncompresed tabular file
DATA <- read.table("data/dmel_mrna.lengc.tbl", header=FALSE);
# otherwise you can run this command
ZZ <- gzfile('data/dmel_mrna.lengc.tbl.gz');
DATA <- read.table(ZZ, header=FALSE);

# just checking the data structure
head(DATA,4);
#         V1   V2    V3
# 1 DQ327735 1883 53.48
# 2 DQ327736  272 55.15
# 3 DQ327737  338 60.36
# 4 DQ327738 1294 50.85
```

We can rename the table columns, so they are more meaningful:

```{.r}
colnames(DATA) <- c("ID","NUClen","GCpct");
head(DATA,4);
#         ID NUClen GCpct
# 1 DQ327735   1883 53.48
# 2 DQ327736    272 55.15
# 3 DQ327737    338 60.36
# 4 DQ327738   1294 50.85
```

Let's calculate some stats on the dataset.

```{.r}
summary(DATA)
#      ID                NUClen            GCpct
# Length:75218       Min.   :   16.0   Min.   :11.43
# Class :character   1st Qu.:  713.2   1st Qu.:47.81
# Mode  :character   Median : 1262.0   Median :51.62
#                    Mean   : 1587.5   Mean   :50.91
#                    3rd Qu.: 2073.0   3rd Qu.:54.82
#                    Max.   :19444.0   Max.   :72.73
```

Now, let's make an histogram:

```{.r}
png(file="images/dmel_hist_length.png");
hist(DATA$NUClen);
dev.off();
png(file="images/dmel_hist_pgc.png");
hist(DATA$GCpct);
dev.off();
```

Or to compare both measures and save into a PNG image...

```{.r}
png(file="images/dmel_plot.png");
plot(DATA$NUClen ~ DATA$GCpct);
dev.off();
```

![Showing GC content versus sequence length](images/dmel_plot.png "Showing GC content versus sequence length")

Just take care that this is not a single-variable distribution but a
two axis scatter-plot, each measure is a continuous random variable
that may fit (or not) a normal distribution. The real data
distribution is shown on the histograms, we can merge then those three
plots... Before doing this, what about getting an improved integration
of images in the report by using a \LaTeX\ floating
environment. Compare result with previous image after compilation, you
can choose to fix that report line.

\input{docs/fig_histograms}


```{.r}
png(file="images/dmel_plot2.png");
def.par <- par();
# preparing a layout grid where to combine different plots
nf <- layout(matrix(c(2,0,1,3), # matrix contents is plot order
                    2, 2, byrow=TRUE),
             c(3,1), c(1,3), # cell relative sizes
             TRUE);
# computing data distribution
xhist <- hist(DATA$GCpct, plot=FALSE);
yhist <- hist(DATA$NUClen, plot=FALSE);
datamax <- max(xhist$counts,
               xhist$counts);
#
# drawing the main plot
par(mar=c(5,5,1,1)); 
plot(DATA$NUClen ~ DATA$GCpct,
      main="",
      xlab="GC%",
      ylab="Sequence length (bp)",
      col="green");
lines(lowess(DATA$NUClen ~ DATA$GCpct),
      col="red", lty="dotted");
mtext(paste("n=",nrow(DATA),sep=""),
      side=3, line=-1);
# drawing x-axis histogram
par(mar=c(0,5,1,1));
barplot(xhist$counts, ylim=c(0, datamax),
        axes=FALSE, space=0, 
        col="steelblue", main="D.melanogaster mRNAs")
# drawing y-axis histogram
par(mar=c(5,0,1,1));
barplot(yhist$counts, xlim=c(0, datamax),
        axes=FALSE, horiz=TRUE, space=0, 
        col="steelblue", main="")
#
par(def.par); # reseting graphical parameters
dev.off();
```

You can try to recode the `dmel_plot2.png` using `ggplot2` `R` library
if you like, but then you must include the corresponding `MarkDown`
code block.

\input{docs/fig_plotcombo}


## Analyzing differences among chromosomes

Let's check if there are changes in the global distributions of mRNAs
length and GC content by chromosome, by strand, or by both; we will
need to merge the chromosome sequence identifiers and strands to
define the corresponding factor variables in our `R` dataframe. We
need a `bed` file relating each mRNA identifier to the chromosome
location, again on the UCSC ftp site from 
["Full database link"](https://hgdownload.soe.ucsc.edu/goldenPath/dm6/database/).

```{.sh}
wget https://hgdownload.soe.ucsc.edu/goldenPath/dm6/database/all_mrna.txt.gz \
     -O $DT/dmel_all_mrna.txt.gz

zcat data/dmel_all_mrna.txt.gz| wc
# 114843 2526546 15008388

# It seems that there are more annotated mRNAs on the genomic sequences
# than mRNA sequences (75218). However, a simple command shows us that
# probably some of those sequences were mapped at multiple genomic locations.
zcat data/dmel_all_mrna.txt.gz | gawk '{print $11}' | sort | uniq -c | wc
#  75130  150260  1274138
```

We need a command-line or script (`Perl`, `Python`, `awk`, ..., as you
like), to merge the strand and chromosome columns from the new
donwloaded file. Just remember to include such script on the
appendices (see \pageref{sec:supplfiles}).

```{.sh}
# a command-line?
# selecting chromosome and strand columns 
zcat dmel_all_mrna.txt.gz | awk '{print $15, $10}' > merged_chr_strand.txt

# sort both file so we can join them 
zcat dmel_mrna.lengc.tbl.gz | sort -k1,1 > sorted_dmel_mrna.lengc.tbl
zcat dmel_all_mrna.txt.gz | sort -k11,11 > sorted_dmel_all_mrna.txt

# join by sequences ID
join -1 1 -2 11 -o 1.1,1.2,1.3,2.10,2.15 sorted_dmel_mrna.lengc.tbl sorted_dmel_all_mrna.txt > merged_mrna.txt
```

Now, the easiest way to compare distributions of lengths and/or GC
content by chromosome, strand or both is to use boxplots. Provide the
`R` commands and include the figure on this report.


```{.r}
library(ggplot2)
merged <- read.table("merged_mrna.txt")
colnames(merged) <- c("ID", "Length", "GC", "Strand", "Chr")

library(dplyr)
# The`Chr` column has chromosome names with additional identifiers, we need to fix this
merged <- merged %>%
  mutate(Chr = gsub("_.*", "", Chr)) %>%  # Clean chromosome labels
  mutate(Chr = as.factor(Chr))            # Convert to factor

# Boxplot of lenghts in all chromosome 
png(file="images/boxplot_leng_chr.png")
ggplot(data = merged, aes(x = Chr, y = log10(Length), fill = Chr)) + 
  geom_boxplot(alpha = 0.5) + 
  geom_jitter(aes(col = Chr), alpha = 0.03, width = 0.25) + 
  labs(x = "All chromosomes", y = "Logarithm of the gene size (bp)") + 
  theme_minimal()
dev.off()

# Boxplot of the GC content in all chromosomes 
png(file="images/boxplot_GC_chr.png")
ggplot(data=merged, aes(x=Chr, y=GC, fill=Chr)) + 
  geom_boxplot(alpha = 0.5) +
  geom_jitter(aes(col = Chr), alpha = 0.03, width = 0.25) +
  labs(x = "All chromosomes", y = "GC content (%)") +
  theme_minimal()
dev.off()

# Lenghts in both strands 
png(file="images/boxplot_leng_strand.png")
ggplot(data=merged, aes(x=Strand, y=log(Length), fill=Strand)) + 
  geom_boxplot(alpha = 0.5) +
  geom_jitter(aes(col = Strand), alpha = 0.03, width = 0.25) +
  labs(x = "DNA Strand", y = "Logarithm of the gene size (bp)")
dev.off() 

# GC content in both strands
png(file="images/boxplot_GC_strand.png")
ggplot(data=merged, aes(x=Strand, y=GC, fill=Strand)) + 
  geom_boxplot(alpha = 0.5) +
  geom_jitter(aes(col = Strand), alpha = 0.03, width = 0.25) +
  labs(x = "DNA Strand", y = "GC content (%)")
dev.off() 

# Lengths of both strands in all chromosomes 
png(file="images/boxplot_leng_chr_strand.png")
ggplot(data=merged, aes(x=Chr, y=log(Length), fill=Strand)) + 
  geom_boxplot(alpha = 0.5) +
  labs(x = "All chromosomes", y = "Logarithm of the gene size (bp)")
dev.off()  

# GC content of both strands in all chromosomes 
png(file="images/boxplot_GC_chr_strand.png")
ggplot(data=merged, aes(x=Chr, y=GC, fill=Strand)) + 
  geom_boxplot(alpha = 0.5) +
  labs(x = "All chromosomes", y = "GC content (%)") 
dev.off()
```


## Further analyses

You can provide here the `bash` and `R` commands to reuse the already
described protocol to analyze another species transcriptome, for
instance the honey bee \Amel\ (["Genome sequence files and select
annotations (2bit, GTF, GC-content, etc)"](https://hgdownload.soe.ucsc.edu/goldenPath/apiMel2/bigZips/)
and ["Full database"](https://hgdownload.soe.ucsc.edu/goldenPath/apiMel2/database/)
links). Then, you can compare the results with respect the fly
transcriptome on the discussion.

```{.sh}
# Your shell commands here

# honey bee files at UCSC genomes repository
wget https://hgdownload.soe.ucsc.edu/goldenPath/apiMel2/bigZips/mrna.fa.gz \
     -O $DT/apimel_mrna.fa.gz
     
zcat $DT/apimel_mrna.fa.gz | head -10
# >GU358185 1
# agtagtgcaggaactgctgtactcataaaatggcgaggatattaaaggta
# ggccaatttggccatacttccacgcgaggtatggaagaatatgtaaaaac
# agtagaacaacgtactcatcatatttcagaaggttctataaaattatgga
# aatctatcactttttttgtcgcttttcctataattggacttgcaatggct
# aattgttatttgaaacatcaagaagaacattcaaaaccgcccccagaatt
# tgttcattatccttatttaaaaatcatgaataagccttttccatggggag
# atggtaaacataca
# >GU358186 1
# gcttttagtttcgcccgggacattagaaccgacgttctaacgtgcaacaa

zcat $DT/apimel_mrna.fa.gz | egrep -c '^>'
# 120637

zcat $DT/apimel_mrna.fa.gz | \
  infoseq -sequence fasta::stdin \
          -noheading -only -name -length -pgc | \
    head -5
# Display basic information about sequences
# GU358185       314    35.35  
# GU358186       719    57.30  
# GU358187       598    52.17  
# GU358188       1117   48.70  
# GU358189       1018   45.58 

zcat $DT/apimel_mrna.fa.gz | \
  infoseq -sequence fasta::stdin \
          -outfile $DT/apimel_mrna.lengc.tbl \
          -noheading -only -name -length -pgc
          
zcat $DT/apimel_mrna.fa.gz | \
  infoseq -sequence fasta::stdin \
          -noheading -only -name -length -pgc | \
    gzip -vc9 - > $DT/apimel_mrna.lengc.tbl.gz
```

```{.sh}
wget https://hgdownload.soe.ucsc.edu/goldenPath/apiMel2/database/all_mrna.txt.gz \
     -O $DT/apimel_all_mrna.txt.gz

zcat apimel_all_mrna.txt.gz| wc
# 127159 2797498 24018653

zcat apimel_all_mrna.txt.gz| gawk '{print $11}' | sort | uniq -c | wc
# 118911  237822 2021443
```

```{.r}
# Your R commands here
DATA <- read.table("apimel_mrna.lengc.tbl", header=FALSE);
head(DATA, 4);
#        V1   V2    V3
# 1 GU358185  314 35.35
# 2 GU358186  719 57.30
# 3 GU358187  598 52.17
# 4 GU358188 1117 48.70
```

```{.r}
colnames(DATA) <- c("ID","NUClen","GCpct");
head(DATA,4);
#         ID NUClen GCpct
# 1 GU358185    314 35.35
# 2 GU358186    719 57.30
# 3 GU358187    598 52.17
# 4 GU358188   1117 48.70
```

```{.r}
summary(DATA)
#      ID                NUClen          GCpct      
# Length:120637      Min.   :   21   Min.   : 0.00  
# Class :character   1st Qu.:  513   1st Qu.:29.80  
# Mode  :character   Median :  937   Median :35.07  
#                    Mean   : 1455   Mean   :35.94  
#                    3rd Qu.: 1938   3rd Qu.:41.87  
#                    Max.   :17525   Max.   :76.36
```

Histogram:
```{.r}
png(file="images/apimel_hist_length.png");
hist(DATA$NUClen);
dev.off();
png(file="images/apimel_hist_pgc.png");
hist(DATA$GCpct);
dev.off();
```
```{.r}
png(file="images/apimel_plot.png");
plot(DATA$NUClen ~ DATA$GCpct);
dev.off();
```

```{.r}
png(file="images/apimel_plot2.png");
def.par <- par();
# preparing a layout grid where to combine different plots
nf <- layout(matrix(c(2,0,1,3), # matrix contents is plot order
                    2, 2, byrow=TRUE),
             c(3,1), c(1,3), # cell relative sizes
             TRUE);
# computing data distribution
xhist <- hist(DATA$GCpct, plot=FALSE);
yhist <- hist(DATA$NUClen, plot=FALSE);
datamax <- max(xhist$counts,
               xhist$counts);
#
# drawing the main plot
par(mar=c(5,5,1,1)); 
plot(DATA$NUClen ~ DATA$GCpct,
      main="",
      xlab="GC%",
      ylab="Sequence length (bp)",
      col="green");
lines(lowess(DATA$NUClen ~ DATA$GCpct),
      col="red", lty="dotted");
mtext(paste("n=",nrow(DATA),sep=""),
      side=3, line=-1);
# drawing x-axis histogram
par(mar=c(0,5,1,1));
barplot(xhist$counts, ylim=c(0, datamax),
        axes=FALSE, space=0, 
        col="steelblue", main="Apis mellifera mRNAs")
# drawing y-axis histogram
par(mar=c(5,0,1,1));
barplot(yhist$counts, xlim=c(0, datamax),
        axes=FALSE, horiz=TRUE, space=0, 
        col="steelblue", main="")
#
par(def.par); # reseting graphical parameters
dev.off();
```


```{.sh}
# selecting chromosome and strand columns 
zcat apimel_all_mrna.txt.gz | awk '{print $15, $10}' > merged_apimel_chr_strand.txt

# sort both file so we can join them 
zcat apimel_mrna.lengc.tbl.gz | sort -k1,1 > sorted_apimel_mrna.lengc.tbl
zcat apimel_all_mrna.txt.gz | sort -k11,11 > sorted_apimel_all_mrna.txt

# join by sequences ID
join -1 1 -2 11 -o 1.1,1.2,1.3,2.10,2.15 sorted_apimel_mrna.lengc.tbl sorted_apimel_all_mrna.txt > merged_apimel_mrna.txt
```

Now, the easiest way to compare distributions of lengths and/or GC
content by chromosome, strand or both is to use boxplots. Provide the
`R` commands and include the figure on this report.


```{.r}
library(ggplot2)
merged <- read.table("merged_apimel_mrna.txt")
colnames(merged) <- c("ID", "Length", "GC", "Strand", "Chr")

library(dplyr)
# The`Chr` column has chromosome names with additional identifiers, we need to fix this
merged <- merged %>%
  mutate(Chr = gsub("_.*", "", Chr)) %>%  # Clean chromosome labels
  mutate(Chr = as.factor(Chr))            # Convert to factor

# Boxplot of lenghts in all chromosome 
png(file="images/apimel_boxplot_leng_chr.png")
ggplot(data = merged, aes(x = Chr, y = log10(Length), fill = Chr)) + 
  geom_boxplot(alpha = 0.5) + 
  geom_jitter(aes(col = Chr), alpha = 0.03, width = 0.25) + 
  labs(x = "All chromosomes", y = "Logarithm of the gene size (bp)") + 
  theme_minimal()
dev.off()

# Boxplot of the GC content in all chromosomes 
png(file="images/apimel_boxplot_GC_chr.png")
ggplot(data=merged, aes(x=Chr, y=GC, fill=Chr)) + 
  geom_boxplot(alpha = 0.5) +
  geom_jitter(aes(col = Chr), alpha = 0.03, width = 0.25) +
  labs(x = "All chromosomes", y = "GC content (%)") +
  theme_minimal()
dev.off()

# Lenghts in both strands 
png(file="images/apimel_boxplot_leng_strand.png")
ggplot(data=merged, aes(x=Strand, y=log(Length), fill=Strand)) + 
  geom_boxplot(alpha = 0.5) +
  geom_jitter(aes(col = Strand), alpha = 0.03, width = 0.25) +
  labs(x = "DNA Strand", y = "Logarithm of the gene size (bp)")
dev.off() 

# GC content in both strands
png(file="images/apimel_boxplot_GC_strand.png")
ggplot(data=merged, aes(x=Strand, y=GC, fill=Strand)) + 
  geom_boxplot(alpha = 0.5) +
  geom_jitter(aes(col = Strand), alpha = 0.03, width = 0.25) +
  labs(x = "DNA Strand", y = "GC content (%)")
dev.off() 

# Lengths of both strands in all chromosomes 
png(file="images/apimel_boxplot_leng_chr_strand.png")
ggplot(data=merged, aes(x=Chr, y=log(Length), fill=Strand)) + 
  geom_boxplot(alpha = 0.5) +
  labs(x = "All chromosomes", y = "Logarithm of the gene size (bp)")
dev.off()  

# GC content of both strands in all chromosomes 
png(file="images/apimel_boxplot_GC_chr_strand.png")
ggplot(data=merged, aes(x=Chr, y=GC, fill=Strand)) + 
  geom_boxplot(alpha = 0.5) +
  labs(x = "All chromosomes", y = "GC content (%)") 
dev.off()
```

![Boxplot of lenghts in all chromosome](images/apimel_boxplot_leng_chr.png "Comparison of all chromosomes based on length")
![Boxplot of the GC content in all chromosomes](images/apimel_boxplot_GC_chr.png "Comparison of all chromosomes based on GC content")

![Lenghts in both strands ](images/apimel_boxplot_leng_strand.png "Comparison of strands based on lengths")
![GC content in both strands](images/apimel_boxplot_GC_strand.png "Comparison of strands based on GC content ")

![Lengths of both strands in all chromosomes](images/apimel_boxplot_leng_chr_strand.png "Comparison of all chromosomes and strands based on lengths")
![GC content of both strands in all chromosomes ](images/apimel_boxplot_GC_chr_strand.png "Comparison of all chromosomes and strands based on GC content")
 

# Discussion
\label{sec:discussion}

__IMPORTANT__ Discuss your results here (around 300 words). And
remember to include in the Appendices section (see page
\pageref{sec:appendices}), any extra script you wrote from this
exercise `bin` folder using the `loadfile` macro. We can take
advantage of the \LaTeX\ referencing capabilities, for instance:
Figure \ref{fig:histograms} on page \pageref{fig:histograms} shows the
GC and sequence lengths histograms distribution for the set of \dmel\
mRNA sequences.


\clearpage

# Appendices
\label{sec:appendices}

## Software

We have used the following versions:

```{.sh}
uname -a
# Linux aleph 5.15.0-117-generic #127-Ubuntu SMP
# Fri Jul 5 20:13:28 UTC 2024 x86_64 x86_64 x86_64 GNU/Linux

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
```


## Supplementary files
\label{sec:supplfiles}


### `conda` environment dependencies for the exercise

\loadfile{environment.yml}{environment.yml}{prg:environmentYML}


### Project specific scripts

\loadfile{an\_script\_example.pl}{bin/an_script_example.pl}{prg:scriptexamplePERL}


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
