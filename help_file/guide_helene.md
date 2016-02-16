@@ -1,241 +0,0 @@
### Read in BAM files and initialize CrisprSet object

Given a set of BAM files with the amplicon sequences of interest mapped to the reference genome, we need to collect a few additional pieces of information
about the guide sequence and define the area around the guide that we want to
summarize the mutation spectrum over.

If you already know the coordinates, these can be typed in or put in a BED file
that can be read in using the **rtracklayer** package.  The import() below command
returns a GRanges object.

```{r, message=FALSE}
library(rtracklayer)
# Represent the guide as a GenomicRanges::GRanges object
gd_fname <- system.file(package="CrispRVariants", "extdata/bed/guide.bed")
gd <- rtracklayer::import(gd_fname)
gd
```

Below, we'll extend the guide region by 5 bases on each side when counting variants.
The guide designed for *ptena* (including PAM) is 23bp and is located on chromosome
`r as.character(seqnames(gd)[1])` from `r start(gd)[1]`-`r end(gd)[1]`.  Note that
the expected cut site (used later for labeling variants), after extension, is 
at base 22 with respect to the start of the guide sequence.

```{r, message=FALSE}
gdl <- resize(gd, width(gd) + 10, fix = "center")
```

With the Bioconductor **BSgenome** packages, the reference sequence itself can be
retrieved directly into a DNAStringSet object.  For other genomes, the reference
sequence can be retrieved from a genome by first indexing the genome with samtools
faidx and then fetching the required region.  Here we are using the GRCHz10
zebrafish genome.  The reference sequence was fetched and saved as follows:

```{r, eval=FALSE}
system("samtools faidx GRCHz10.fa.gz")

reference=system(sprintf("samtools faidx GRCHz10.fa.gz %s:%s-%s", 
                         seqnames(gdl)[1], start(gdl)[1], end(gdl)[1]) , 
                 intern = TRUE)[[2]]

# The guide is on the negative strand, so the reference needs to be reverse complemented
reference=Biostrings::reverseComplement(Biostrings::DNAString(reference))
save(reference, file = "ptena_GRCHz10_ref.rda")
```

We'll load the previously saved reference sequence.

```{r}
ref_fname <- system.file(package="CrispRVariants", "extdata/ptena_GRCHz10_ref.rda")
load(ref_fname)
reference
```

Note the NGG sequence (here, TGG) is present with the 5 extra bases on the end.

To allow easy matching to experimental condition (e.g., useful for colour labeling)
and for subsetting to experiments of interest, we often organize the list of BAM
files together with accompanying metadata in a machine-readable table beforehand:

```{r, message=FALSE}
# The metadata and bam files for this experiment are included with CrispRVariants
library("gdata")
md_fname <- system.file(package="CrispRVariants", "extdata/metadata/metadata.xls")
md <- gdata::read.xls(md_fname, 1)
md

# Get the bam filenames from the metadata table
bam_dir <- system.file(package="CrispRVariants", "extdata/bam")
bam_fnames <- file.path(bam_dir, md$bamfile)

# check that all files exist
all( file.exists(bam_fnames) )
```

## Creating a CrisprSet 
The next step is to create a CrisprSet object, which is the container that stores
the relevant sequence information, alignments, observed variants and their frequencies.

```{r, message=FALSE}
# Note that the zero point (target.loc parameter) is 22
crispr_set <- readsToTarget(bam_fnames, target = gdl, reference = reference,
                            names = md$Short.name, target.loc = 22)
crispr_set

# The counts table can be accessed with the "variantCounts" function
vc <- variantCounts(crispr_set)
print(class(vc))
```

You can see that in the table of variant counts, variants are summarised by the
location of their insertions and deletions with respect to the target site.
Non-variant sequences and sequences with a single nucleotide variant (SNV) but
no insertion or deletion (indel) are displayed first, followed by the indel
variants from most to least frequent  For example, the most frequent
non-wild-type variant, "-1:4D" is a 4 base pair deletion starting 1 base
upstream of the zero point.

## Creating summary plots of variants

We want to plot the variant frequencies along with the location of the guide
sequence relative to the known transcripts.  If you do this repeatedly for the
same organism, it is worthwhile to save the database in a local file and read
in with loadDb(), since this is quicker than retrieving it from UCSC
(or Ensembl) each time.

We start by creating a transcript database of Ensembl genes.  The gtf was
downloaded from Ensembl version 81.  We first took a subset of just the genes
on chromosome 17 and then generated a transcript database.

```{r, engine='bash', eval=FALSE}
# Extract genes on chromosome 17 (command line)
# Note that the Ensembl gtf does not include the "chr" prefix, so we add it here 
gtf=Danio_rerio.GRCz10.81.gtf.gz
zcat ${gtf} | awk '($1 == 17){print "chr"$0}' > Danio_rerio.GRCz10.81_chr17.gtf
```

```{r, eval = FALSE}
# In R
library(GenomicFeatures)
gtf_fname <- "Danio_rerio.GRCz10.81_chr17.gtf"
txdb <- GenomicFeatures::makeTxDbFromGFF(gtf_fname, format = "gtf")
saveDb(txdb, file= "GRCz10_81_chr17_txdb.sqlite")
```
We now load the the previously saved database 

```{r, echo=FALSE, message=FALSE}
library(GenomicFeatures)
txdb_fname <- system.file("extdata/GRCz10_81_ptena_txdb.sqlite", 
                          package="CrispRVariants")
txdb <- loadDb(txdb_fname)
```

plotVariants() is a wrapper function that groups together a plot of the
transcripts of the gene/s overlapping the guide (optional),
CrispRVariants::plotAlignments(), which displays the alignments of the
consensus variant sequences to the reference, and
CrispRVariants::plotFreqHeatmap(), which produces a table of the variant
counts per sample, coloured by either their counts or percentage contribution
to the variants observed for a given sample.  If a transcript database is supplied,
the transcript plot is  annotated with the guide location. Arguments for
plotAlignments() and plotFreqHeatmap() can be passed to plotVariants() as lists
named plotAlignments.args and plotFreqHeatmap.args, respectively.

```{r}
# The gridExtra package is required to specify the legend.key.height 
# as a "unit" object.  It is not needed to call plotVariants() with defaults
library(gridExtra)

# Match the clutch id to the column names of the variants
group <- md$Group
```


```{r, fig.width = 8.5, fig.height = 9.5, message = FALSE, fig.cap = "(Top) schematic of gene structure showing guide location (left) consensus sequences for variants (right) variant counts in each embryo."}
p <- plotVariants(crispr_set, txdb = txdb, gene.text.size = 8, 
    row.ht.ratio = c(1,8), col.wdth.ratio = c(4,2),
    plotAlignments.args = list(line.weight = 0.5, ins.size = 2, 
                               legend.symbol.size = 4),
    plotFreqHeatmap.args = list(plot.text.size = 3, x.size = 8, group = group, 
                                legend.text.size = 8, 
                                legend.key.height = grid::unit(0.5, "lines"))) 
```

The plotVariants() options set the text size of the transcript plot annotation
(gene.text.size) and the relative heights (row.ht.ratio) and widths
(col.wdth.ratio) of the plots. 

The plotAlignments arguments set the symbol size in the figure (ins.size)
and in the legend (legend.symbol), the line thickness for the (optional)
annotation of the guide region and cleavage site (line.weight). 

For plotFreqHeatmap we define an grouping variable for colouring the x-axis
labels (group), the size of the text within the plot (plot.text.size) and on
the x-axis (x.size) and set the size of the legend text (legend.text.size).

## Calculating the mutation efficiency

The mutation efficiency is the number of reads that include an insertion or
deletion.  Chimeric reads and reads containing single nucleotide variants near
the cut site may be counted as variant reads, non-variant reads, or excluded
entirely.  See the help page for the function **mutationEfficiency** for more
details.

We can see in the plot above that the control sample includes a variant
sequence 6:1D, also present in sample ptena 4.  We will exclude all sequences
with this variant from the efficiency calculation.  We also demonstrate below
how to exclude particular variants. 


```{r}
# Calculate the mutation efficiency, excluding indels that occur in the "control" sample
# and further excluding the "control" sample from the efficiency calculation
eff <- mutationEfficiency(crispr_set, filter.cols = "control", exclude.cols = "control")
eff

# Suppose we just wanted to filter particular variants, not an entire sample.
# This can be done using the "filter.vars" argument
eff2 <- mutationEfficiency(crispr_set, filter.vars = "6:1D", exclude.cols = "control")

# The results are the same in this case as only one variant was filtered from the control
identical(eff,eff2)
```

We see above that sample ptena 4 has an efficiency of 80%, i.e. 4 variant sequences,
plus one sequence "6:1D" which is counted as a non-variant sequence as it also
occurs in the control sample.

## Plot chimeric alignments

When deciding whether chimeric alignments should be considered as variant sequences,
it can be useful to plot the frequent chimeras.  


```{r}
ch <- getChimeras(crispr_set, sample = "ptena 4")

# Confirm that all chimeric alignments are part of the same read
length(unique(names(ch))) == 1

# Set up points to annotate on the plot
annotations <- c(resize(gd, 1, fix = "start"), resize(gd, 1, fix = "end"))
annotations$name <- c("ptena_start", "ptena_end")

plotChimeras(ch, annotations = annotations)

```

Here we see the read aligns as two tandem copies of the region chr17:23648420-23648656.
The endpoint of each copy is not near the guide sequence.  We do not consider this a
genuine mutation, so we'll recalculate the mutation efficiency excluding the chimeric
reads and the control variant as before.

```{r}
mutationEfficiency(crispr_set, filter.cols = "control", exclude.cols = "control",
                   include.chimeras = FALSE)
```

We see that the mutation effiency for "ptena 4" is now 75%, i.e. 3 genuine variant
sequences, 1 sequence counted as "non-variant" because it occurs in the control,
and the chimeric read excluded completely.
