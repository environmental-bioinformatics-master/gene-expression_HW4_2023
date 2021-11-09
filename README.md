# 2021_HW4B_Gene-Expression-2

## DUE: 11:59 PM, 15 November

From the first part of this homework, you should have a tab-separated .tsv table with a different contig in each row and a different sample in each column, giving the raw number of reads that uniquely mapped to each contig for each sample.

Next step is to determine which, if any, of these contigs are differentially expressed between your two treatment groups. The most popular DE analysis approaches are based in R, which we are not covering extensively in the class. These programs have extensive step-by-step help and vignettes available online. If you're not conversant in R, please make good use of these resources, and come talk to the instructors for additional help if needed.

Copy this document and change the name to hw3b_answers_\[LASTNAME\].md, and reply in the document as prompted. We'll also ask you to create and save an R script and some results files to be pushed with your HW repo. Everything we want you to submit is listed at the bottom of this document.

Please use the differential expression R package assigned here:

DESeq2: Walter, Annaliese, Sabrina, Matthew, Katy, Caroline
edgeR: Evie, Mira, Kate, Jane, Max, Sun

Note: If you're not formally enrolled but want to do this homework, please pick whichever DE package strikes your fancy.

You are welcome to interface with R however you prefer for this homework: on Poseidon (via Jupyter if desired) or on your personal computer (via RStudio if desired). Likewise, feel free to work with packages either through conda or directly through R.

Reminder: An extended note on making figures in R on a cluster: Do NOT try to make figures directly using R on Poseidon. By default, R prints to screen…and the HPC doesn’t have a screen. (It will not print to your local computer, since it’s not running locally.) If you want to make a figure, first open an "empty" pdf file, execute the command(s) that would typically print a figure to screen, and then close the pdf like this:

pdf(file = FILENAME.pdf)
EXECUTE_COMMANDS_TO_MAKE_FIGURE
dev.off()

To view the resulting pdf, either copy it to your local computer with scp, or keep an Jupyter notebook open in your working directory and view the file in that notebook.

## Part 1: Background and setup

How are you running R? (Poseidon, personal computer, etc.)
> Answer:

How did you install your DE packages? (conda or install.packages)
> Answer:

What DE approach are you using?
> Answer:

What resources did you use to figure out how to run it (please give URLs, etc)?
> Answer:

Is it recommended for any particular type of data / analysis?
> Answer:

After installing your DE package(s) and loading them via the library(PACKAGE) command, type sessionInfo(). This will show you everything that's been loaded, and its version number. Other attached packages gives details on everything loaded directly (rather than called as a dependency only). What packages are listed under "other attached packages"? Include package name and version number.

> Answer:

## Step 2: Load data and prep for DE analysis

To test for differential expression of eye genes at night vs during the day, you'll need to associate samples with their capture metadata. I've provided a file called `HW4_sample_metadata.txt` listing this information for each sample. Note that this files lists the SRR numbers only, while the sample names in your count file may be named slightly differently. These names will need to be identical between the two files for downstream processing. Change the sample names in your counts table to match those in the sample treatment file. You can use any approach you like. How did you do this?

>Answer:

For your homework, create a new R script called `hw4b_de_\[LASTNAME\].R`. In it, put all of the commands you use to analyse your data, starting with loading R packages. We don't need to be able to run this script from the command line, but we should be able to copy your commands into R and replicate your analysis.

First, you need to load your data in a format that your DE program can read. Most of the programs we're using require data in matrix format. This should load your count data in an appropriate format for downstream work:
```
counts = as.matrix(read.csv("COUNT_TABLE_FILE", sep="\t", row.names="Contig"))
```

Check your matrix with dim(counts) and head(counts), to make sure it loaded properly. What is the result of your dim(counts) command?
>Answer:

Now, load the treatment data:
treatment = read.csv("HW3_sample_treatments.txt", sep = "\t", row.names=1, header = T)

...and check that it has loaded properly.

Your samples are probably not in the same order in the counts and treatment files, but they should be. Re-order the columns in your counts table to match the (numerical) ordering in the treatment file like this:
counts = counts[, rownames(treatment)]

Check that this re-ordering worked with head(counts).
Step 3: Begin your DE analysis

For some approaches, you may want to pre-filter your data set to remove genes with very low expression. This slims down the data by removing genes whose expression is too low to give you any useful information on differential expression, which is helpful when you adjust your p-values to account for multiple tests. Various approaches handle low-expression filtering differently, and some explicitly do not recommend it; use the approach recommended for your DE package.

Did you filter or pre-filter your data to remove low-expressed transcripts? If so, what cutoff did you use and why? If not, why not?
Answer:

What command(s) did you use to filter your data?

How many transcripts are left in your data set after filtering (or not)?
Answer:

Following the guidance for your DE package, load your data into your DE program and associate your treatment information with it. What commands did you use to do this?

Note: Depending on your approach, you may need to load the data before filtering it.

The first step in some approaches is to normalize your data. This compensates for differences in sequencing depth across samples, sequencing biases, etc. This can be done differently depending on the DE analysis program you’re using. Some programs have this type of normalization baked-in, while you have to do it explicitly in others.

What, if any, normalization did you do on your data?
Answer:

What metric did you use to normalize your data? What does this mean?
Answer:

What command(s) did you use to normalize your data?

Step 3: Conduct DE analysis

Now it’s time to look for differential expression between your two groups. Follow the guidance in your program’s help files on this. Your output should be a table of log-fold changes and p-values for your control:treatment comparison for each gene. (It may include other columns as well.) Make sure p-values are adjusted to alpha = 0.05.

What commands did you use to test differential expression in your data?

How does your DE program model distribution of read counts (e.g. Poisson, negative binomial, etc)?
Answer:

How many different transcripts are you using in this analysis? (Count just those with enough reads mapping to be informative the final DE analysis, excluding any removed from earlier filtering or determined by the program to have low counts.)
Answer:

Write your DE table to a tab-separated file like this:
write.table(RESULTS_TABLE, file = "hw3_de-table_[LASTNAME].tsv", sep = "\t")

Push your results .tsv file as part of your homework.
Step 4: Explore your results

How many genes are differentially expressed at an adjusted standard p-value or equivalent (p < 0.05)?
Answer:

How many genes are differentially expressed at an adjusted stringent p-value or equivalent (p < 0.01)?
Answer:

How many genes show a log-fold change > 2?
Answer:

How many genes show a log-fold change > 5?
Answer:

What are the 20 most-differentiated transcripts based on p-value?
Answer:

What are the 10 transcripts most upregulated in the treatment vs. control, based on log-fold change?
Answer:

What are the 10 transcripts most downregulated in the treatment vs. control, based on log-fold change?
Answer:

Note: All programs have a lot of neat built-in plotting and other functions to explore your data. Feel free to check them out - they are often very helpful - but we're not asking you to include any plots in this homework.

About how long did this homework take you?
Answer:
For your homework, please push to GitHub:

    hw3b_answers_[LASTNAME].md: An annotated copy of this readme file including your answers.
    hw3_de_[LASTNAME].R: A file containing all of the R commands you used in your analysis.
    hw3_de-table_[LASTNAME].tsv: A tab-separated file of your DE results.
