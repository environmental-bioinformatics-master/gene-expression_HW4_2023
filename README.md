# 2021_HW4B_Gene-Expression-2

## DUE: 11:59 PM, 15 November

From the first part of this homework, you should have a tab-separated .tsv table with a different contig in each row and a different sample in each column, giving the raw number of reads that uniquely mapped to each contig for each sample.

The next step is to determine which, if any, of these contigs are differentially expressed between your two treatment groups. The most popular DE analysis approaches are based in R, which we are not covering extensively in the class. These programs have extensive step-by-step help and vignettes available online. If you're not conversant in R, please make good use of these resources, and come talk to the instructors for additional help if needed.

Copy this document and change the name to hw3b_answers_\[LASTNAME\].md, and reply in the document as prompted. We'll also ask you to create and save a jupyter notebook in R, a modified R script, and some results files and figures to be pushed with your HW repo. Everything we want you to submit is listed at the bottom of this document.

Please use the differential expression R package assigned here:

**DESeq2**: Walter, Annaliese, Sabrina, Matthew, Katy, Caroline

**edgeR**: Evie, Mira, Kate, Jane, Max, Sun

Note: If you're not formally enrolled but want to do this homework, please pick whichever DE package strikes your fancy.

For this part of the homework, modify the r_jupyter yaml file to create a new environment called `HW4B_env`. Make sure to run the extra line of code to set up an R kernel in Jupyter: `R -e 'IRkernel::installspec()'`. Install all of the R programs you use in the environment directly in R, and then spin up a new jupyter notebook in R called `HW4B_diffex.ipynb`. Run all your R commands for differential expression analysis in this notebook, and push it as part of the HW.

Reminder: An extended note on making figures in R on a cluster: Do NOT try to make figures directly using R on Poseidon. By default, R prints to screen…and the HPC doesn’t have a screen. (It will not print to your local computer, since it’s not running locally.) If you want to make a figure, first open an "empty" pdf file, execute the command(s) that would typically print a figure to screen, and then close the pdf like this:

```
pdf(file = FILENAME.pdf)
EXECUTE_COMMANDS_TO_MAKE_FIGURE
dev.off()
```

To view the resulting pdf, either copy it to your local computer with scp, or keep an Jupyter notebook open in your working directory and view the file in that notebook.

## Part 1: Background and setup

What DE approach are you using?
> Answer:

Is it recommended for any particular type of data / analysis?
> Answer:

What resources did you use to figure out how to run it (please give URLs, etc)?
> Answer:

What commands did you use to install R packages in your environment?

```
```

After installing your DE package(s) and loading them in your new notebook via the library(PACKAGE) command, type sessionInfo(). This will show you everything that's been loaded, and its version number. Other attached packages gives details on everything loaded directly (rather than called as a dependency only). What packages are listed under "other attached packages"? Include package name and version number.
> Answer:

## Step 2: Load data and prep for DE analysis

To test for differential expression of eye genes at night vs during the day, you'll need to associate samples with the depth at whih they were captured (your experimental "treatments"). I've provided a file called `hw4_sample_metadata.csv` listing all the SRA information given for each sample. Note that this files lists the SRR numbers only, while the sample names in your count file may be named slightly differently. These names **have to** be identical between the two files for downstream processing. If necessary, change the sample names in your counts table to match those in the sample metadata file. You can use any approach you like. How did you do this?

>Answer:

First, you need to load your data in a format that your DE program can read. Most of the programs we're using require data in matrix format. This should load your count data in an appropriate format for downstream work:

```
counts = as.matrix(read.csv("COUNT_TABLE_FILE", sep="\t", row.names="Contig"))
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

For some approaches, you may want to pre-filter your data set to remove genes with very low expression. This slims down the data by removing genes whose expression is too low to give you any useful information on differential expression, which is helpful when you adjust your p-values to account for multiple tests. Various approaches handle low-expression filtering differently, and some explicitly do not recommend it; use the approach recommended for your DE package. Note: Depending on your approach, you may need to load the data into the analysis program before filtering it.

Did you filter or pre-filter your data to remove low-expressed transcripts? If so, what cutoff did you use and why? If not, why not?
>Answer:

What command(s) did you use to filter your data?

How many transcripts are left in your data set after filtering (or not)?
>Answer:

Following the guidance for your DE package, load your data into your DE program and associate your treatment (depth_factor) information with it. What commands did you use to do this?
>Answer:

## Step 3: Conduct DE analysis

Now it’s time to look for differential expression between your two groups. Follow the guidance in your program’s help files on this. Your output should be a table of log-fold changes and p-values for your control:treatment comparison for each gene. (It may include other columns as well.) Make sure p-values are adjusted to alpha = 0.05.

What commands did you use to test differential expression in your data?
>Answer:

How does your DE program model distribution of read counts (e.g. Poisson, negative binomial, etc)?
>Answer:

How many different transcripts are you using in this analysis? (Count just those with enough reads mapping to be informative the final DE analysis, excluding any removed from earlier filtering or determined by the program to have low counts.)
>Answer:

Write your DE table to a comma-separated file like this:

```
write.csv(RESULTS_TABLE, file = "hw4b_de-table_[LASTNAME].csv", quote = FALSE)
```

The `quote = FALSE` is important - without this, R will write the contig names in quotation marks, which won't match the quotation-mark-free contig names in the GO annotation file, and the GO_MWU program you'll use for functional analysis will refuse to believe that they are the same names. (Isn't compatibility fun?)

Push your results .csv file as part of your homework.

## Step 4: Explore your results

How many genes are differentially expressed at an adjusted standard p-value or equivalent (p < 0.05)?
>Answer:

How many genes are differentially expressed at an adjusted stringent p-value or equivalent (p < 0.01)?
>Answer:

How many genes show a log-fold change > 2?
>Answer:

>How many genes show a log-fold change > 5?
Answer:

What are the 10 most-differentiated transcripts based on p-value?
>Answer:

What are the 5 transcripts most upregulated in the treatment vs. control, based on log-fold change?
>Answer:

What are the 5 transcripts most downregulated in the treatment vs. control, based on log-fold change?
>Answer:

Why aren't the *same* genes identified using p-values and log2-fold change values?
>Answer:

Note: All programs have a lot of neat built-in plotting and other functions to explore your data. Feel free to run some if you're keen - they are often very helpful - but you don't have to include any differential expression plots in this homework.


## Step 5: Functional enrichment

Now that you have a list of p-values and log2fold-change data, what does it all mean? Let's explore functional enrichment a bit, using an R-based approach called `GO_MWU`: https://github.com/z0on/GO_MWU

This program expects a graphical output and it's a pain to work around this, so copy any necessary files to your local computer and work with GO_MWU in R on your local system. (You may also need to install perl, which is required by the program.)

Clone the GO_MWU repo and follow the directions for running this analysis. The developers helpfully include a "template" file, `GO_MWU.R`, for running in R. Copy this file as `hw4b_go-mwu_\[LASTNAME\].R` and modify this file to apply to your data.

Run GO_MWU on your data using adjusted p-values as your "continuous measure of change". Run this analysis 3 times: once for each GO category (Molecular Function, Cellular Component, Biological Process). There are some parameters you can "tune" to help highlight the most important functional categories without overwhelming the output with a zillion lines - "largest", "smallest", and "clusterCutHeight". Read about these in the program documentation, and play around with them to find a set of parameters that you think represents the data accurately and clearly.

What parameters did you pick for each GO analysis? Why?
>Molecular Function:

>Biological Process:

>Cellular Compartment:


Save the resulting figures as `hw4b_go-mwu-res_\[GO-CATEGORY\]_\[LASTNAME\].pdf`. Scale the figure if needed to ensure it's legible when opened.

How do you interpret these graphs? Why might this approach (ranked list) be better than enrichment based on a p-value cutoff for this analysis?

>Answer:

Using your biological knowledge and all the analyses you've run on these data, what kind of physiological changes do you think are important to diel vertical migration in these shrimp? Please explain briefly (a paragraph, two max):

>Answer:

What caveats do you have about these analyses? If you had the opportunity to conduct one new experiment to explore this question in more detail, what would you propose (a paragraph, two max)?

>Answer:


About how long did this homework take you?
>Answer:


For your homework, please push to GitHub:

- `hw4b_answers_[LASTNAME].md`: An annotated copy of this readme file including your answers.
- `HW4B_diffex.ipynb`: A jupyter notebook containing all of the R commands you used in your differential expression analysis
- `hw4b_de-table_[LASTNAME].csv`: The results file from your differential expression analysis (table with contig name, p-value, log2-fold change, etc).
- `hw4b_go-mwu_\[LASTNAME\].R`: GO_MWU running file, modified for your analysis (any one of the GO categories is fine here, no need to include 3 copies).
- `hw4b_go-mwu-res_MF_\[LASTNAME\].pdf`: Graph of Molecular Function enrichment analysis.
- `hw4b_go-mwu-res_BP_\[LASTNAME\].pdf`: Graph of Biological Process enrichment analysis.
- `hw4b_go-mwu-res_CC_\[LASTNAME\].pdf`: Graph of Cellular Compartment enrichment analysis.
   
