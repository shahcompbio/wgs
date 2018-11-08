library(TitanCNA)

version <- '0.1.5'

args <- commandArgs(TRUE)

id <- args[1]
tc_het_file <- args[2]
cnfile <- args[3]
map <- args[4]
numClusters <- as.numeric(args[5])
numCores <- as.numeric(args[6])
ploidy <- as.numeric(args[7])
outfile <- args[8]
outparam <- args[9]
myskew <- as.numeric(args[10])
boolEstPloidy <- args[11]
n_zero <- as.numeric(args[12])
normEstMeth <- args[13]
maxI <- as.numeric(args[14])
pseudo_counts =  as.numeric(args[15])
txn_exp_len = as.numeric(args[16])
txn_z_strength = as.numeric(args[17])
alphaK <- as.numeric(args[18])        #prior for events; default: 15000
alphaHigh <- as.numeric(args[19])     #prior for extreme events; default: 15000
maxCN <- as.numeric(args[20])         #maximum number of copies to use
sym <- args[21]
outobj <- args[22]
genometype <- args[23]
chrom <- args[24]
yThreshold <- as.numeric(args[25])
maxDep <- as.numeric(args[26])
chrom <- eval(parse(text=chrom)) 

message('Running TITAN...')

#### LOAD DATA ####
data <- loadAlleleCounts(tc_het_file, symmetric=sym, genomeStyle=genometype)

#### LOAD PARAMETERS ####
message('titan: Loading default parameters')
params <- loadDefaultParameters(copyNumber=maxCN, numberClonalClusters=numClusters, skew=myskew, symmetric=sym)
params$ploidyParams$phi_0 <- ploidy
params$normalParams$n_0 <- n_zero
    
# #### GC AND MAPPABILITY CORRECTION ####
message('titan: Reading GC content and mappability corrected read counts ...')
cnData <- read.delim(cnfile, header=TRUE, stringsAsFactors=FALSE, sep="\t")

#### READ COPY NUMBER FROM HMMCOPY FILE ####
message('titan: Extracting read depth...')
logR <- getPositionOverlap(data$chr, data$posn,cnData)

data$logR <- log(2^logR)
rm(logR, cnData)

#### FILTER DATA FOR DEPTH, MAPPABILITY, NA, etc ####
mScore <- as.data.frame(wigToRangedData(map))
mScore <- getPositionOverlap(data$chr, data$posn, mScore[,-4])

#### Check if Chromosomes Have been provided

if (is.null(chrom)) {
chrom <- unique(sort(data$chr))
}

# check if sample is Female or number of datapoints is very small.
if (genometype=='UCSC'){
    if (NROW(filterData(data, c('chrY'), minDepth=10, maxDepth=maxDep, map=mScore, mapThres=0.8)) > yThreshold){
        data <- filterData(data, chrom, minDepth=10, maxDepth=maxDep, map=mScore, mapThres=0.8)
    } else {
        data <- filterData(data, chrom[which(chrom!='chrY')], minDepth=10, maxDepth=maxDep, map=mScore, mapThres=0.8)
    }
} else{
    if (NROW(filterData(data, c('Y'), minDepth=10, maxDepth=maxDep, map=mScore, mapThres=0.8)) > yThreshold){
        data <- filterData(data, chrom, minDepth=10, maxDepth=maxDep, map=mScore, mapThres=0.8)
    } else {
        data <- filterData(data,chrom[which(chrom!='Y')], minDepth=10, maxDepth=maxDep, map=mScore, mapThres=0.8)
    }
}

#### MODEL SELECTION USING EM (FWD-BACK) TO SELECT NUMBER OF CLUSTERS ####
library(doMC)
registerDoMC(cores=numCores)

##### RUN USING EM ALGORITHM ######
K <- length(params$genotypeParams$rt)
params$genotypeParams$alphaKHyper <- rep(alphaK,K)
if (sym) { highStates <- c(1,7:K) } else { highStates <- c(1,11:K) }
params$genotypeParams$alphaKHyper[highStates] <- alphaHigh
convergeParams <- runEMclonalCN(data, gParams=params$genotypeParams,
								nParams=params$normalParams,
                                pParams=params$ploidyParams, sParams=params$cellPrevParams,
                                maxiter=maxI, maxiterUpdate=1500, txnExpLen=txn_exp_len,
                                txnZstrength=txn_z_strength,                                                
                                useOutlierState=FALSE, normalEstimateMethod=normEstMeth,
                                estimateS=TRUE, estimatePloidy=boolEstPloidy,
                                pseudoCounts=pseudo_counts)
    
#### COMPUTE OPTIMAL STATE PATH USING VITERBI ####
#options(cores=1)
optimalPath <- viterbiClonalCN(data, convergeParams)

#### PRINT RESULTS TO FILES ####
tryCatch({
    results <- outputTitanResults(data, convergeParams, optimalPath, filename=outfile,posteriorProbs=FALSE, subcloneProfiles=TRUE)
    outputModelParameters(convergeParams, results, outparam)

    convergeParams$ploidy = ploidy
    convergeParams$numClusters = numClusters
 
    save(convergeParams, results, file=paste(outobj))
},
error = function(err){
    print('setting subcloneprofiles to False and retrying due to error:')
    print(err)
    results <- outputTitanResults(data, convergeParams, optimalPath, filename='outfile', posteriorProbs=FALSE, subcloneProfiles=FALSE)
    outputModelParameters(convergeParams, results, outparam)
    save(convergeParams, results, file=paste(outobj))
})
