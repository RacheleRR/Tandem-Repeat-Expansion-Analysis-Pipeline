# DBSCAN 

pkgs <- c("dbscan", "data.table", "parallel", "doSNOW")

for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, repos = "https://cloud.r-project.org")
  }
}

# Load necessary libraries
library(dbscan)
library(data.table)
library(parallel)
library(doSNOW)

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
outpath <- args[1]
infile <- args[2]
samplelist <- args[3]
coreNumber <- as.integer(args[4])


# Use the parsed arguments in the script
params <- list(
  outpath = outpath,
  infile = infile,
  samplelist = samplelist
)

# Debugging: Print input file
print(paste("Debug: params$infile is set to:", params$infile))

# Set working directory
# setwd(params$outpath)
print(paste("Current working directory is:", getwd()))

# Check if the input file exists
if (!file.exists(params$infile)) {
  stop(paste("Error: The file does not exist:", params$infile))
}

# Debugging: Displaying paths
cat("DEBUG: The input file passed to R is:", params$infile, "\n")
cat("DEBUG: The output directory is:", params$outpath, "\n")

set.seed(2000)

outpath <- params$outpath
if (length(grep("\\/$", outpath)) == 0) {
  outpath <- paste0(outpath, "/")
}

if (!dir.exists(outpath)) {
  dir.create(outpath, recursive = TRUE)
}

infile <- params$infile
message(sprintf("reading %s##", infile))
dt <- fread(infile)
dt <- data.frame(dt)
dt$varid <- paste(dt$V1, dt$V2, dt$V3, sep="_")

# Get count per sample
sample.count <- paste(dt$V7, collapse = ",")
sample.count <- gsub(":", ",", sample.count)
sample.count <- strsplit(sample.count, ",")[[1]]
sample.count <- matrix(sample.count, ncol = 2, byrow = TRUE)
sample.count <- data.frame(table(sample.count[, 1]), stringsAsFactors = FALSE)

if ("samplelist" %in% names(params)) {
  samples <- readLines(params$samplelist)
  sample.count <- sample.count[sample.count$Var1 %in% samples, ]
} else {
  samples <- sample.count$Var1
}

sd <- sd(sample.count$Freq)
mean <- mean(sample.count$Freq)
outlier.count <- sum(sample.count$Freq > mean + 3 * sd | sample.count$Freq < mean - 3 * sd)

if (outlier.count > 0) {
  message(sprintf("Remove %s outliers", outlier.count))
  samples <- samples[!samples %in% sample.count$Var[which(sample.count$Freq > mean + 3 * sd | sample.count$Freq < mean - 3 * sd)]]
}

message(sprintf("Samples included in the analysis are in %sclean.sample.txt", outpath))
writeLines(as.character(sample.count$Var1), sprintf("%sclean.sample.txt", outpath))

norm.dist <- round(rnorm(length(samples), mean = 1, sd = 0.25), digits = 1)
norm.dist <- sort(norm.dist)

time.start <- Sys.time()

cl <- makeCluster(coreNumber-1)
registerDoSNOW(cl)

dt.out <- foreach (i = 1:nrow(dt), .combine=rbind) %dopar% {
  tmp <- gsub(":|,", "#", dt[i, 7])
  tmp <- unlist(strsplit(tmp, "#")[[1]])
  tmp <- data.frame(matrix(tmp, ncol = 2, byrow = T), stringsAsFactors = F)
  tmp[, 2] <- as.numeric(tmp[, 2])
  tmp <- tmp[tmp$X1 %in% samples, ]
  if(nrow(tmp) == 0){
    dt.here <- data.frame("repeatID" = "", "motif" = "", "outliers" = "", "size" = "", "ref" = "",
                          "chr" = "", "start" = "", "end" = "")[-1, ]
  }else{
    if(nrow(tmp) < length(norm.dist)){
      dt.tmp <- norm.dist[-c(1:nrow(tmp))]
      sampleIDs <- c(paste0("Sim", 1:length(dt.tmp)), tmp$X1)
      dt.tmp <- c(dt.tmp, tmp$X2)
    }else{
      sampleIDs <- tmp$X1
      dt.tmp <- tmp$X2
    }
    
    ref <- as.numeric(names(which.max(table(dt.tmp))))
    range <-  max(2*ref, quantile(dt.tmp, 0.95, na.rm = T) - quantile(dt.tmp, 0.05, na.rm = T))
    
    scan <- dbscan::dbscan(matrix(dt.tmp), eps = range, minPts = ceiling(log2(length(samples))))
  
    if(length(unique(scan$cluster)) == 1 | sum(scan$cluster == 0) == 0){
      cutoff <- Inf
    }else{
      cutoff <- max(dt.tmp[scan$cluster != 0])
      cutoff <- ifelse(cutoff < 2, 2, cutoff)
    }
    
    if(sum(dt.tmp > cutoff) > 0){
      sim.samples <- grep("Sim", sampleIDs)
      if(length(sim.samples) > 0){
        dt.tmp <- dt.tmp[-sim.samples]
        sampleIDs <- sampleIDs[-sim.samples]
      }
      
      outliers <- paste(sampleIDs[dt.tmp > cutoff], collapse=";")
      size <- paste(dt.tmp[dt.tmp > cutoff], collapse=";")
      dt.here <- data.frame("repeatID" = dt$varid[i], "motif" = dt$V4[i], outliers, size, ref,
                            "chr" = dt$V1[i], "start" = dt$V2[i], "end" = dt$V3[i])
    }else{
      dt.here <- data.frame("repeatID" = "", "motif" = "", "outliers" = "", "size" = "", "ref" = "",
                            "chr" = "", "start" = "", "end" = "")[-1, ]
    }
  }

  dt.here
}

stopCluster(cl) 

output_file <- file.path(outpath, "DBSCAN_expansions.tsv")

write.table(
  dt.out,
  output_file,
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)

message(paste("DBSCAN results written to:", output_file))
