### Convert AB1-format Sanger sequences to FASTQ

If the raw data comes as
AB1 (Sanger) format.  To convert AB1 files to FASTQ, we use ab1ToFastq(), which is a wrapper for functions in **sangerseqR** package with additional quality score trimming.

[See User guide on github](https://github.com/markrobinsonuzh/CrispRVariants/blob/master/vignettes/user_guide.Rmd)

### Map the FASTQ reads

We use FASTQ format because it is the major format used by most genome alignment algorithms.