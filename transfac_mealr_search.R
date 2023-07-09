#!/usr/bin/Rscript
#
# Copyright (c) 2022 geneXplain GmbH, Wolfenbüttel, Germany
#
# Author: Philip Stegmaier, philip.stegmaier@genexplain.com
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY 
# OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT
# LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
# FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO
# EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE 
# LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, 
# WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
# ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
# OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#

library(geneXplainR)
library(optparse)
library(stringr)

get_cmd_options <- function() {
    parser = OptionParser()

    parser <- add_option(parser, c("-i", "--fasta"),
                         type     = "character",
                         default  = NULL,
                         help     = "Input FASTA file")
    
    parser <- add_option(parser, c("-n", "--name"),
                         type     = "character",
                         default  = "mealr_search_input",
                         help     = "Designated name of sequence in the platform workspace. This is the name of the sequence that will be uploaded into the output folder. [default= %default]")
    
    parser <- add_option(parser, c("-g", "--gene"),
                         type    = "character",
                         default = "",
                         help    = "Input gene symbol")
    
    parser <- add_option(parser, c("-d", "--species"),
                         type    = "character",
                         default = "Human (Homo sapiens)",
                         help    = "Species to associate with specified gene symbol [default= %default]. Strings that can be specified for the species are listed in the 'Convert table' tool's 'Species' option (https://platform.genexplain.com/#de=analyses/Methods/Data%20manipulation/Convert%20table)")
    
    parser <- add_option(parser, c("-m", "--genome"),
                         type    = "character",
                         default = "Ensembl 104.38 Human (hg38)",
                         help    = "Platform name of reference genome for specified gene [default= %default]. Available genomes are listed in MATCH(TM) tool's sequence sources (https://platform.genexplain.com/#de=analyses/Methods/Site%20analysis/TRANSFAC(R)%20Match(TM)%20for%20tracks)")
    
    parser <- add_option(parser, c("-a", "--start"),
                         type    = "integer",
                         default = -1000,
                         help    = "TSS-relative sequence start, for specified gene symbol [default= %default]")
    
    parser <- add_option(parser, c("-b", "--end"),
                         type    = "integer",
                         default = -100,
                         help    = "TSS-relative sequence end, for specified gene symbol [default= %default]")
    
    parser <- add_option(parser, c("-r", "--result"),
                         type     = "character",
                         default  = "mealr_search_result",
                         help     = "Designated name of the search result folder in the platform workspace. This folder is created within the output folder (-p/--platform_out) [default= %default]")
    
    parser <- add_option(parser, c("-o", "--out"),
                         type    = "character",
                         default = "transfac_mealr_search_result",
                         help    = "Local output folder [default= %default]")
    
    parser <- add_option(parser, c("-y", "--prob"),
                         type    = "double",
                         default = 0.99,
                         help    = "Filter results for predictions with specified probability cutoff [default= %default]")
    
    parser <- add_option(parser, c("-f", "--factor"),
                         type    = "character",
                         default = "",
                         help    = "Filter results for specified factor gene symbol. Multiple can be specified separated by a pipe symbol (|)")
    
    parser <- add_option(parser, c("-c", "--cell"),
                         type    = "character",
                         default = "",
                         help    = "Filter results for specified cell source. Multiple can be specified separated by a pipe symbol (|)")
    
    parser <- add_option(parser, c("-t", "--tissue"),
                         type    = "character",
                         default = "",
                         help    = "Filter results for specified tissue source. Multiple can be specified separated by a pipe symbol (|)")
    
    parser <- add_option(parser, c("-e", "--profile"),
                         type    = "character",
                         default = "databases/TRANSFAC(R) 2022.2/Data/profiles/vertebrate_human_p0.001",
                         help    = "TRANSFAC(R) reference profile. Score thresholds from this profile are used for binding site prediction. This profile should cover the vertebrate set of PWMs")
    
    parser <- add_option(parser, c("-p", "--platform_out"),
                         type    = "character",
                         default = "",
                         help    = "Output folder in the platform workspace. A string like 'data/Projects/<project name>/Data/<subfolder(s)>'. The folder is created if it does not exist. The parent folder is expected to be present")
    
    parser <- add_option(parser, c("-s", "--server"),
                         type    = "character",
                         default = "https://platform.genexplain.com",
                         help    = "The platform server URL to connect with [default= %default]")
    
    parser <- add_option(parser, c("-u", "--user"),
                         type    = "character",
                         default = "",
                         help    = "Platform user for login")
    
    parser <- add_option(parser, c("-q", "--password"),
                         type    = "character",
                         default = "",
                         help    = "Platform user's password for login")
    
    parser <- add_option(parser, c("-w", "--passfile"),
                         type    = "character",
                         default = "",
                         help    = "Source file with platform user's password for login to avoid giving the password on the commandline. The source contains a variable 'gxpass' to which the password is assigned like 'gxpass <- \"password\"'")
    
    opt = parse_args(parser)
    
    if (nchar(opt$fasta) == 0 && nchar(opt$gene) == 0) {
        stop("Error: missing input FASTA file or gene (one of options -i/--fasta or -g/--gene required)")
    }
    
    if (nchar(opt$platform_out) == 0) {
        stop("Error: missing platform output folder (option -p/--platform_out)")
    }
    
    if (nchar(opt$user) == 0) {
        stop("Error: missing platform username (option -u/--user)")
    }
    
    if (nchar(opt$password) == 0 && nchar(opt$passfile) == 0) {
        stop("Error: missing platform password (one of options -q/--password or -w/--passfile required)")
    }
    
    opt
}


#
# Creates specified folder in platform folder it it doesn't exist
#
prepare_platform_folder <- function(platformOut) {
    pout <- str_split(platformOut, fixed("/"))[[1]]
    parentFolder <- paste(pout[1:length(pout)-1], collapse = "/")
    if (!gx.isElement(parentFolder, pout[length(pout)])) {
        gx.createFolder(parentFolder, pout[length(pout)])
    }
}


split_multiple_args <- function(arg) {
    str_trim(str_split(arg, fixed("|"))[[1]], "both")
}


#
# Parses FASTA file and extracts first sequence
#
get_sequence <- function(fasta) {
    con   <- file(fasta, open="r")
    lines <- readLines(con)
    close(con)
    inSeq = FALSE
    sid <- ""
    seq <- ""
    for (line in lines){
        if (str_sub(line, 1, 1) == ">") {
            if (inSeq) { break }
            sid <- str_sub(line, 2)
            inSeq <- TRUE
        } else if (inSeq) {
            seq <- paste0(seq, str_replace(line, fixed(" "), ""))
        }
    }
    c(sid, seq)
}


#
# Creates for specified sequence
#
write_temp_fasta <- function(seq, fasta) {
    lines <- c(paste0(">", seq[1]), seq[2])
    outfile <- file(fasta, 'w')
    writeLines(lines, outfile)
    close(outfile)
}


#
# Extracts table rows that contain terms
# in specified column
#
filter_result_table <- function(table, terms, col) {
    ftable <- table
    terms <- terms[nchar(terms) > 0]
    if (length(terms) > 0 && nchar(terms[1]) > 0) {
        ftable <- table[which(table[, col] %in% terms),]
    }
    ftable
}


#
# Extracts value of specified attribute field from a GTF row
#
get_gtf_field <- function(g, name, isNumeric = FALSE) {
    anno <- strsplit(strsplit(g, "; ", fixed = TRUE)[[1]], ' "', fixed = TRUE)
    v <- str_replace_all(anno[[which(sapply(anno, function(x){ x[1] == name}))]][2], '"', "")
    if (isNumeric) {
        v <- as.numeric(v)
    }
    v
}



cmdArgs <- get_cmd_options()

cmdArgs$factor <- split_multiple_args(cmdArgs$factor)
cmdArgs$cell   <- split_multiple_args(cmdArgs$cell)
cmdArgs$tissue <- split_multiple_args(cmdArgs$tissue)

#
# Create local output folder
#
dir.create(cmdArgs$out, showWarnings = FALSE)

#
# Sign into platform workspace
#
if (nchar(cmdArgs$passfile) > 0) {
    source(cmdArgs$passfile)
    cmdArgs$password <- gxpass
}
gx.login(cmdArgs$server, cmdArgs$user, cmdArgs$password)

#
# Create folder in platform workspace if it doesn't exist
#
prepare_platform_folder(cmdArgs$platform_out)

inputName <- ""
inputPath <- ""

if (nchar(cmdArgs$gene) > 0) {
    #
    # Check sequence coordinates
    #    
    pstart <- cmdArgs$start
    if (pstart == 0) {
        pstart = 1
    }
    pend <- cmdArgs$end
    if (pend == 0) {
        pend = 1
    }
    if (pstart > pend) {
        stop("Error: sequence start position (-a/--start) > sequence end position (-b/--end)")
    }
    
    slen <- pend + 1 - pstart
    if (pstart < 0 && pend > 0) {
        slen <- pend - pstart
    }
    if (slen < 50) {
        stop("Error: sequence too short. Sequence length >= 50 required")
    }
    
    #
    # Import specified gene
    #
    inputFile <- paste0(cmdArgs$gene, ".txt")
    inputFilepath <- file.path(cmdArgs$out, inputFile)
    write.table(data.frame(Symbol = c(cmdArgs$gene)), file = inputFilepath, quote = FALSE, row.names = FALSE)
    gx.importTable(inputFilepath,cmdArgs$platform_out,
                   cmdArgs$gene, columnForID = "Symbol", tableType = "Genes: Gene symbol", species = cmdArgs$species)
    
    # Short break after data import
    Sys.sleep(1)
    
    #
    # Convert gene symbol to Ensembl gene
    #
    inputGenePath    <- paste0(cmdArgs$platform_out, "/", cmdArgs$gene)
    inputEnsemblPath <- paste0(cmdArgs$platform_out, "/", cmdArgs$gene, " Ensembl")
    gx.analysis("Convert table",
                list(sourceTable = inputGenePath,
                     species     = cmdArgs$species,
                     sourceType  = "Genes: Gene symbol",
                     targetType  = "Genes: Ensembl",
                     outputTable = inputEnsemblPath),
                TRUE, TRUE)
    
    #
    # Create track of genomic region around TSS (promoter)
    #
    inputPath <- paste0(cmdArgs$platform_out, "/", cmdArgs$gene, " sequence")
    gx.analysis("Gene set to track",
                list(sourcePath = inputEnsemblPath,
                     species    = cmdArgs$species,
                     from       = pstart,
                     to         = pend,
                     destPath   = inputPath),
                TRUE, TRUE)
    inputName <- cmdArgs$gene
} else {
    #
    # Get (first) sequence from FASTA file
    #
    seq <- get_sequence(cmdArgs$fasta)
    
    if (nchar(seq[1]) < 50) {
        stop("Error: sequence too short. Sequence length >= 50 required")
    }
    
    #
    # Makes single sequence FASTA file for upload
    #
    fastaFile <- paste0(cmdArgs$name, ".fasta")
    tempFasta <- file.path(cmdArgs$out, fastaFile)
    write_temp_fasta(seq, tempFasta)

    #
    # Import input sequence
    #
    print("Importing sequence into platform workspace")
    gx.import(tempFasta, cmdArgs$platform_out, "Fasta format (*.fasta)")
    
    inputName <- cmdArgs$name
    inputPath  <- paste0(cmdArgs$platform_out, "/", inputName)
}

# Short break after data import
Sys.sleep(1)

#
# Prepare platform result folder path
#
searchName <- cmdArgs$result
outputPath <- paste0(cmdArgs$platform_out, "/", searchName)

#
# Scan sequence with TRANSFAC(R) MEALR models
#
print("Starting combinatorial regulation analysis")
gx.analysis("TRANSFAC(R) MEALR combinatorial regulation analysis",
            params <- list(sequencePath = inputPath,
                           asClassifier = FALSE,
                           stepSize     = 50,
                           scanMode     = "Best hit",
                           accuracyCut  = 0.85,
                           output       = outputPath),
            TRUE, TRUE)

#
# Export MEALR result table
#
resultTable <- paste0(outputPath, "/MEALR search result")
localResultTable <- file.path(cmdArgs$out, "MEALR_search_result.tsv")
gx.export(resultTable, target.file = localResultTable)

#
# Filter result table for hit probability, factor, tissue, cell source
#
mealrResult <- read.table(localResultTable, header = TRUE, sep = "\t", check.names = FALSE)

filteredResult <- subset(mealrResult, Probability >= cmdArgs$prob)

filteredResult <- filter_result_table(filteredResult, cmdArgs$factor, "Factor gene")
filteredResult <- filter_result_table(filteredResult, cmdArgs$tissue, "Tissue source")
filteredResult <- filter_result_table(filteredResult, cmdArgs$cell, "Cell source")

#
# Store and import filtered result table
#
localFilteredResultTable <- file.path(cmdArgs$out, "MEALR_search_result_filtered.tsv")

write.table(filteredResult, file = localFilteredResultTable, quote = FALSE, sep = "\t", row.names = FALSE)

filteredResultName <- "MEALR search result filtered"
filteredResultPath <- paste0(outputPath, "/MEALR search result filtered")

gx.importTable(localFilteredResultTable, outputPath, filteredResultName,
               delim = "Tab", headerRow = 1, dataRow = 2, tableType = "Unspecified")

# Short break after data import
Sys.sleep(1)

#
# Export MEALR result sequence track to GTF file
#
resultTrack <- paste0(outputPath, "/MEALR search result track")
localResultTrack <- file.path(cmdArgs$out, "MEALR_search_result.gtf")
gx.export(resultTrack, exporter = "Gene Transfer Format (*.gtf)", exporter.params = list(includeHeader = FALSE), target.file = localResultTrack)

#
# Filter result track for hit probability, factor, tissue, cell source
#
trackResult <- read.table(localResultTrack, header = FALSE, sep = "\t", quote = "")
trackResult$Probability <- as.vector(sapply(trackResult[,9], function(x){ get_gtf_field(x, "probability", TRUE) }))
trackResult <- subset(trackResult, Probability >= cmdArgs$prob)
trackResult$FactorGene <- as.vector(sapply(trackResult[,9], function(x){ get_gtf_field(x, "factorGene", FALSE) }))
trackResult <- filter_result_table(trackResult, cmdArgs$factor, "FactorGene")
trackResult$TissueSource <- as.vector(sapply(trackResult[,9], function(x){ get_gtf_field(x, "tissueSource", FALSE) }))
trackResult <- filter_result_table(trackResult, cmdArgs$tissue, "TissueSource")
trackResult$CellSource <- as.vector(sapply(trackResult[,9], function(x){ get_gtf_field(x, "cellSource", FALSE) }))
trackResult <- filter_result_table(trackResult, cmdArgs$cell, "CellSource")

#
# Store and import filtered result track
#
localFilteredTrackTable <- file.path(cmdArgs$out, "MEALR_search_result_track_filtered.gtf")

write.table(trackResult[,1:9], file = localFilteredTrackTable, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

gx.import(localFilteredTrackTable, outputPath, "Gene Transfer Format (*.gtf)", list(dbSelector = cmdArgs$genome))

filteredTrackPath <- paste0(outputPath, "/MEALR_search_result_track_filtered")

# Short break after data import
Sys.sleep(1)

#
# Extract PWMs from selected MEALR models
#
filteredPwmsPath <- paste0(outputPath, "/MEALR search result PWMs")

gx.analysis("Extract TRANSFAC(R) PWMs from combinatorial regulation analysis",
            params <- list(mealrOutputPath = filteredResultPath,
                           output = filteredPwmsPath),
            TRUE, TRUE)

#
# Create MATCH(TM) profile with extracted PWMs using cutoffs
# from reference profile
#
mealrProfilePath <- paste0(outputPath, "/MEALR search result PWMs profile")

gx.analysis("Create profile from site model table",
            list(table = filteredPwmsPath,
                 profile = cmdArgs$profile,
                 outputProfile = mealrProfilePath),
            TRUE, TRUE)

#
# Predict binding sites with selected PWMs within locations
# of selected models
#
matchPath <- paste0(outputPath, "/Filtered MEALR PWM binding sites")

gx.analysis("TRANSFAC(R) Match(TM) for tracks", 
            list(sequencePath = filteredTrackPath,
                 profilePath  = mealrProfilePath,
                 withoutDuplicates = TRUE,
                 ignoreCore = TRUE,
                 output = matchPath),
            TRUE, TRUE)

#
# Download and store Match(TM) result in GTF file
#
localMatchTrack <- file.path(cmdArgs$out, "MEALR_PWM_binding_sites.gtf")
gx.export(matchPath, exporter = "Gene Transfer Format (*.gtf)", exporter.params = list(includeHeader = FALSE), target.file = localMatchTrack)

gx.logout()
