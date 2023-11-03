# 2023_HW4_Gene-Expression

## DUE: 11:59 PM, 14 November

With your Snakemake homework, you should have generated a tab-separated .tsv table with a different contig in each row and a different sample in each column, giving the raw number of reads that uniquely mapped to each contig for each sample.

The next step is to determine which, if any, of these contigs are differentially expressed between your two treatment groups. The most popular DE analysis approaches are based in R, which we are not covering extensively in the class. These programs have extensive step-by-step help and vignettes available online. If you're not conversant in R, please make good use of these resources, and come talk to the instructors for additional help if needed.

Copy this document and change the name to hw4_de_answers_\[LASTNAME\].md, and reply in the document as prompted. We'll also ask you to create and save an R script file, a modified R script, and some results files and figures to be pushed with your HW repo. Everything we want you to submit is listed at the bottom of this document.

You will be running most of your analyses in R. You have options! You can work on your local computer (in RStudio if you prefer), or on the HPC directly, or on the HPC in a Jupyter notebook with an R kernel. Whatever you choose, I want you to make a well-commented file containing *all of your R commands* for differential expression - this can be in *.R (R script) or *.ipynb (Jupyter notebook) format.

Reminder: An extended note on making figures in R on a cluster: Do NOT try to make figures directly using R on Poseidon. By default, R prints to screen…and the HPC doesn’t have a screen. (It will not print to your local computer, since it’s not running locally.) If you want to make a figure, first open an "empty" pdf file, execute the command(s) that would typically print a figure to screen, and then close the pdf like this:

```
pdf(file = FILENAME.pdf)
EXECUTE_COMMANDS_TO_MAKE_FIGURE
dev.off()
```

To view the resulting pdf, either copy it to your local computer with scp, or keep an Jupyter notebook open in your working directory and view the file in that notebook.

## Part 1: Background and setup

You should already have DESeq2 installed on the HPC. If you don't, or if you want to work locally, install it. After loading DESeq2 (and any other necessary packages) into R via the library(PACKAGE) command, type sessionInfo(). This will show you everything that's been loaded, and its version number. Other attached packages gives details on everything loaded directly (rather than called as a dependency only). What packages are listed under "other attached packages"? Include package name and version number.
> Answer:

What resources did you use to figure out how to run DESeq2 (please give URLs, etc)?
> Answer:


## Step 2: Load data and prep for DE analysis

To test for differential expression of eye genes at night vs during the day, you'll need to associate samples with the depth at which they were captured (your experimental "treatments"). I've provided a file called `hw4_sample_metadata.csv` listing all the SRA information given for each sample. Note that this files lists the SRR numbers only, while the sample names in your count file may be named slightly differently. These names **have to** be identical between the two files for downstream processing. If necessary, change the sample names in your counts table to match those in the sample metadata file. You can use any approach you like (does not need to be in R). Did you need to do this? If so, how did you do this?

>Answer:

First, you need to load your data in a format that DESeq2 can recognize - like many DE programs, it likes matrix format. This basic command should load your count data in an appropriate format for downstream work:

```
counts = as.matrix(read.csv("COUNT_TABLE_FILE", sep="\t", row.names="CONTIG_COL_NAME"))
```

Check your matrix with dim(counts) and head(counts), to make sure it loaded properly. What is the result of your dim(counts) command?
>Answer:

Now, load the treatment data:

```
treatment = read.csv("hw4_sample_metadata.csv", row.names=1, header = T)
```

...and check that it has loaded properly. The factor we're going to compare is `depth_factor`, where I've defined samples as "deep" or "shallow". (R really does not appreciate spaces and dashes as in the original `depth` column.)

Your samples are probably not in the same order in the counts and treatment files, but they should be. Re-order the columns in your counts table to match the (numerical) ordering in the treatment file like this:

```
counts = counts[, rownames(treatment)]
```

Check that this re-ordering worked with head(counts).

## Step 3: Begin your DE analysis

For some approaches, you may want to pre-filter your data set to remove genes with very low expression. This slims down the data by removing genes whose expression is too low to give you any useful information on differential expression, which is helpful when you adjust your p-values to account for multiple tests. Various approaches handle low-expression filtering differently, and some explicitly do not recommend it; always check the approach recommended for your DE package.

Load your data into DESeq2 and associate your treatment (depth_factor) information with it. (This may require some additional prep of your data files.)

How many transcripts are in your raw data set?
>Answer:

Filter your data to remove all transcripts with < 20 reads mapped to them.

How many transcripts are left in your data set after filtering?
>Answer:

## Step 3: Conduct DE analysis

Now it’s time to look for differential expression between your two groups. Follow the DESeq2 guidance on this. Your output should be a table of log-fold changes and p-values for your control:treatment comparison for each gene. (It may include other columns as well.) Make sure p-values are adjusted to alpha = 0.05.

How does DESeq2 model distribution of read counts (e.g. Poisson, negative binomial, etc)?
>Answer:

How many different transcripts are you using in this analysis? (Count just those with enough reads mapping to be informative the final DE analysis, excluding any removed from earlier filtering or determined by the program to have low counts.)
>Answer:

Write your DE table to a comma-separated file like this:

```
write.csv(RESULTS_TABLE, file = "hw4_de-table_[LASTNAME].csv", quote = FALSE)
```

The `quote = FALSE` is important - without this, R will write the contig names in quotation marks, which won't match the quotation-mark-free contig names in the GO annotation file, and the GO_MWU program you'll use for functional analysis will refuse to believe that they are the same names. (Isn't compatibility fun?)

Push your results .csv file as part of your homework.

## Step 4: Explore your results

Remember to update your R  / jupyter script with any commands you used to answer these questions!

How many genes are differentially expressed at an adjusted standard p-value or equivalent (p < 0.05)?
>Answer:

How many genes are differentially expressed at an adjusted stringent p-value or equivalent (p < 0.01)?
>Answer:

How many genes show a log-fold change > 2 or < -2?
>Answer:

How many genes show a log-fold change > 5 or < -5?
>Answer:

What are the 10 most-differentiated transcripts based on p-value?
>Answer:

What are the 5 transcripts most upregulated in the treatment vs. control, based on log-fold change?
>Answer:

What are the 5 transcripts most downregulated in the treatment vs. control, based on log-fold change?
>Answer:

Why aren't the *same* genes identified using p-values and log2-fold change values?
>Answer:

Note: DESeq2 has a lot of neat built-in plotting and other functions to explore your data. Feel free to run some if you're keen - they are often very helpful - but you don't have to include any differential expression plots in this homework.


## Step 5: Functional enrichment

Now that you have a list of p-values and log2fold-change data, what does it all mean? Let's explore functional enrichment a bit, using an R-based approach called `GO_MWU`: https://github.com/z0on/GO_MWU

This program expects a graphical output and it's a pain to work around this, so copy any necessary files to your local computer and work with GO_MWU in R on your local system. (You may also need to install perl, which is required by the program.)

Clone the GO_MWU repo and follow the directions for running this analysis. The developers helpfully include a "template" file, `GO_MWU.R`, for running in R. Copy this file as `hw4_go-mwu_LASTNAME.R` and modify this file to apply to your data. You will need to download a files in the symlinked directory from HW3 to run GO_MWU.

What reference file did you need from your symlinked directory?
>Answer:

Run GO_MWU on your data using adjusted p-values as your "continuous measure of change". You will first need to modify the data file you produced in the prior step. Use the command line to do this.

What command did you use to modify you DE table to run GO_MWU?
>Answer:

Copy your modified data file (for running GO_MWU) into your homework repo.

Run GO_MWU 3 times: once for each GO category (Molecular Function, Cellular Component, Biological Process). There are some parameters you can "tune" to help highlight the most important functional categories without overwhelming the output with a zillion lines - "largest", "smallest", and "clusterCutHeight". Read about these in the program documentation, and play around with them to find a set of parameters that you think represents the data accurately and clearly.

What parameters did you pick for each GO analysis? Why?
>Molecular Function:

>Biological Process:

>Cellular Compartment:


Save the resulting figures as `hw4_go-mwu-res_GO-CATEGORY_LASTNAME.pdf`. Scale the figure if needed to ensure it's legible when opened.

How do you interpret these graphs? Why might this approach (ranked list) be better than enrichment based on a p-value cutoff for this analysis?

>Answer:

Using your biological knowledge and all the analyses you've run on these data, what kind of physiological changes do you think are important to diel vertical migration in these shrimp? Please explain briefly (a paragraph, two max):

>Answer:

What caveats do you have about these analyses? If you had the opportunity to redesign the sequencing strategy for this project to explore this question in more detail, what would you propose (a paragraph, two max)?. Assume you have access to a good number of samples collected as per the paper, and a reasonable (though not overly generous) sequencing budget. 

>Answer:


About how long did this homework take you?
>Answer:


For your homework, please push to GitHub:

- `hw4_answers_LASTNAME.md`: An annotated copy of this readme file including your answers.
- `hw4_diffex_LASTNAME.R` OR `HW4_diffex_LASTNAME.ipynb`: A commented R script OR jupyter notebook containing all of the R commands you used in your differential expression analysis
- `hw4_de-table_LASTNAME.csv`: The results file from your differential expression analysis (table with contig name, p-value, log2-fold change, etc).
- `hw4_go-mwu_LASTNAME.R`: GO_MWU running file, modified for your analysis (any one of the GO categories is fine here, no need to include 3 copies).
- The modified DE data file you used to run GO_MWU (no required naming convention, as long as it's the same name you use for it in `hw4_go-mwu_LASTNAME.R`).
- `hw4_go-mwu-res_MF_LASTNAME.pdf`: Graph of Molecular Function enrichment analysis.
- `hw4_go-mwu-res_BP_LASTNAME.pdf`: Graph of Biological Process enrichment analysis.
- `hw4_go-mwu-res_CC_LASTNAME.pdf`: Graph of Cellular Compartment enrichment analysis.
   
