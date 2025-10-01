
The code of how I purturb the model is one [github](https://gitlab.com/uniluxembourg/lcsb/BioCore/disease_maps/machine-learning-for-boolean-networks/genetic-algorithms-bingyu-jiang-internship/-/blob/bf8fadcb280a0b5cd17c56d90ca0b9735458ac65/02.Pertub_model.R). 

Following is the code to run the process: 

```r
file="./data/ToyModel/10_Modified/ToyModel.RData"
# file = "data/DREAMmodel/DreamModel.RData"
# data(CNOlistDREAM,package="CellNOptR")
# data(DreamModel,package="CellNOptR")
# file = "data/ToyModel/ToyModel.RData"
load(file)
pknmodel <- mod_model
library("MEIGOR") # will load CellNOptR
data("CellNOptR_example", package="MEIGOR") # same ToyModel as in CellNOptR
cnolist <- CNOlist(cnolist_cellnopt)
model <- preprocessing(data = cnolist, model = pknmodel)
t <- system.time(res <- gaBinaryT1(CNOlist = cnolist,model =  model, verbose=FALSE))
optModel <- cutModel(model, res$bString);
# A small bug of cutModel function, it did not reconvert the notMat matrix.

if (ncol(optModel$notMat) == 0) {{
    n <- nrow(optModel$interMat)
    m <- ncol(optModel$interMat)
    optModel$notMat <- matrix(0, nrow=n, ncol=m)
    rownames(optModel$notMat) <- rownames(optModel$interMat)
    colnames(optModel$notMat) <- colnames(optModel$interMat)
    
    for (j in seq_along(optModel$reacID)) {{
        rid <- optModel$reacID[j]
        print(paste("Processing reaction ID:", rid))
        is_not <- startsWith(rid, "!")
        rid_clean <- sub("^!", "", rid)
        parts <- strsplit(rid_clean, "=")[[1]]
        src <- parts[1]
        tgt <- parts[2]
        if (is_not) {{
            print(paste("Adding NOT reaction:", src, "->", rid))
            optModel$notMat[src, rid] <- 1
        }}
    }}
}}

library(here)
source(here::here("tools", "functions.R"))

sif_fname = "ToyModel_10.sif"
boolnet_fname = "ToyModel_10.bnet"

writeSIF(optModel, file=sif_fname, overwrite=TRUE) # A function of CellNOptR

SIFToBoolNet(sifFile     = sif_fname,
            boolnetFile = boolnet_fname,
            CNOlist     = cnolist,
            model       = optModel,
            fixInputs   = FALSE,
            preprocess  = TRUE,
            ignoreAnds  = TRUE)
```
will yield the following error:
```txt
character(0)
Error in cat(genes[i], ", ", geneStrings[i], "\n", sep = "") : 
  argument 3 (type 'list') cannot be handled by 'cat'
```

The `SIFToBoolNet` is located in [`tools/functions.R`](https://gitlab.com/uniluxembourg/lcsb/BioCore/disease_maps/machine-learning-for-boolean-networks/genetic-algorithms-bingyu-jiang-internship/-/blob/bf8fadcb280a0b5cd17c56d90ca0b9735458ac65/tools/functions.R), which works fine for the original model, but not for the perturbed model. 

```r
library("CellNOptR")
pknmodel <- readSIF("data/ToyModel/0_Modified/ToyModel.sif")
midas_file = "data/ToyModel/0_Modified/ToyModel.csv"
md<- readMIDAS(midas_file, verbose=TRUE)
cnolist <- makeCNOlist(md, subfield=FALSE)
model <- preprocessing(data = cnolist, model = pknmodel)

library(here)
# Source helper functions
source(here::here("tools", "functions.R"))
t <- system.time(resILP <- ilpBinaryT1New(
                    CNOlist(cnolist), model, md,
                    cplexPath="/home/users/bjiang/CPLEX_Studio2211/cplex/bin/x86-64_linux/cplex", 
                           sizeFac = 0.0001, mipGap=0, relGap=0, 
                           timelimit=3600, 
                           method = "quadratic", numSolutions = 100, 
                           limitPop = 500, poolIntensity = 0, 
                           poolReplace = 2))

print("Running ILP optimization...");
cnolist <- CNOlist(cnolist);
resILPAll <- list();
exclusionList <- NULL;
cnolistReal <- cnolist;       
writeMIDAS(CNOlist = cnolist, filename = "tempMD.csv", 
    timeIndices = c(1, 2), overwrite = TRUE);
md <- readMIDAS(MIDASfile = "tempMD.csv");
file.remove("tempMD.csv");
cnolist <- compatCNOlist(object = cnolist); 
options(scipen = 10); 
cnolist <- makeCNOlist(md,FALSE);

t <- system.time(resILP <- CellNOptR:::createAndRunILP(
                    model, md, cnolistReal,
                    cplexPath="~/CPLEX_Studio2211/cplex/bin/x86-64_linux/cplex", accountForModelSize = TRUE, 
                           sizeFac = 0.0001, mipGap=0, relGap=0, 
                           timelimit=3600, 
                           method = "quadratic", numSolutions = 100, 
                           limitPop = 500, poolIntensity = 0, 
                           poolReplace = 2))


```