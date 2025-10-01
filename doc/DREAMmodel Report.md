```r
summarizeDreamData <- function(model, CNOl) {
  # 1) Network (BoolNet) summary
  n_nodes       <- length(model$namesSpecies)
  n_reactions   <- length(model$reacID)
  n_activations <- sum(model$interMat == 1)
  n_inhibitions <- sum(model$notMat   == 1)
  
  cat("=== Network summary ===\n")
  cat("• Nodes (species)         :", n_nodes, "\n")
  cat("• Reactions (edges)       :", n_reactions, "\n")
  cat("  – Activating edges      :", n_activations, "\n")
  cat("  – Inhibitory edges      :", n_inhibitions, "\n\n")
  
  # 2) Experimental data (CNOlist) summary
  n_cues       <- length(CNOl$namesCues)
  n_stimuli    <- length(CNOl$namesStimuli)
  n_inhibitors <- length(CNOl$namesInhibitors)
  n_signals    <- length(CNOl$namesSignals)
  n_times      <- length(CNOl$timeSignals)
  n_samples    <- nrow(CNOl$valueCues)
  
  cat("=== Experimental data summary ===\n")
  cat("• Cues (inputs measured)    :", n_cues, "\n")
  cat("• Stimuli applied           :", n_stimuli, "\n")
  cat("• Inhibitors applied        :", n_inhibitors, "\n")
  cat("• Signals read out          :", n_signals, "\n")
  cat("• Time‐points per experiment:", n_times, "\n")
  cat("• Total experiments (rows)  :", n_samples, "\n")
}
library(CellNOptR)
midas_file <- "data/TCell/TCell.csv"
sif_file <- "data/TCell/TCell.sif"
pknmodel <- readSIF(sif_file)
cnolist <- makeCNOlist(readMIDAS(midas_file), subfield = FALSE)  
summarizeDreamData(pknmodel, cnolist)

data(CNOlistDREAM,package="CellNOptR")
data(DreamModel,package="CellNOptR")

data(CNOlistToy,package="CellNOptR")
data(ToyModel,package="CellNOptR")
summarizeDreamData(DreamModel, CNOlistDREAM)
summarizeDreamData(ToyModel, CNOlistToy)

```
=== Network summary ===
• Nodes (species)         : 40 
• Reactions (edges)       : 58 
  – Activating edges      : 58 
  – Inhibitory edges      : 2 

=== Experimental data summary ===
• Cues (inputs measured)    : 8 
• Stimuli applied           : 4 
• Inhibitors applied        : 4 
• Signals read out          : 7 
• Time‐points per experiment: 2 
• Total experiments (rows)  : 25


> CNOlistDREAM$namesCues
[1] "igf1"  "il1a"  "tgfa"  "tnfa"  "ikk"   "mek12" "pi3k"  "p38"  
> CNOlistDREAM$namesInhibitors
[1] "ikk"   "mek12" "pi3k"  "p38"  
> CNOlistDREAM$namesStimuli
[1] "igf1" "il1a" "tgfa" "tnfa"
> CNOlistDREAM$namesSignals
[1] "akt"   "erk12" "ikb"   "jnk12" "p38"   "hsp27" "mek12"
> DreamModel
$reacID
 [1] "prak=hsp27"    "map3k7=mkk7"   "map3k7=nik"    "map3k7=mkk3"  
 [5] "map3k7=mkk6"   "map3k7=mkk4"   "pip3=pdk1"     "pak=raf1"     
 [9] "rac=map3k1"    "rac=pak"       "raf1=mek12"    "mek12=erk12"  
[13] "irs1=pi3k"     "nik=ikk"       "map3k1=mkk7"   "map3k1=ikk"   
[17] "map3k1=mkk4"   "grb2=sos"      "mkk4=p38"      "mkk4=jnk12"   
[21] "igf1=igfr"     "mkk7=jnk12"    "igfr=irs1"     "igfr=shc"     
[25] "sos=ras"       "ikk=ikb"       "!akt=raf1"     "akt=cot"      
[29] "!akt=pak"      "traf6=sitpec"  "traf6=map3k7"  "il1r=traf6"   
[33] "p38=prak"      "shc=grb2"      "mkk6=p38"      "tgfa=egfr"    
[37] "traf2=ask1"    "traf2=map3k7"  "sitpec=map3k1" "pi3k=rac"     
[41] "pi3k=pip3"     "ras=pi3k"      "ras=raf1"      "ras=map3k1"   
[45] "egfr=shc"      "egfr=grb2"     "egfr=pi3k"     "tnfr=traf2"   
[49] "tnfr=pi3k"     "il1a=il1r"     "cot=ikk"       "mkk3=p38"     
[53] "mkk3=prak"     "erk12=prak"    "pdk1=akt"      "tnfa=tnfr"    
[57] "ask1=mkk4"     "ask1=mkk7"    

$namesSpecies
 [1] "prak"   "map3k7" "pip3"   "pak"    "rac"    "raf1"   "mek12"  "irs1"  
 [9] "nik"    "map3k1" "grb2"   "mkk4"   "igf1"   "mkk7"   "igfr"   "sos"   
[17] "ikk"    "akt"    "traf6"  "il1r"   "p38"    "shc"    "mkk6"   "tgfa"  
[25] "traf2"  "sitpec" "pi3k"   "ras"    "egfr"   "tnfr"   "il1a"   "cot"   
[33] "mkk3"   "erk12"  "pdk1"   "tnfa"   "ask1"   "hsp27"  "jnk12"  "ikb" 

