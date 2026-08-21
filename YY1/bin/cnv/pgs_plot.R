library(karyoploteR)

args <- commandArgs(TRUE)
if (length(args) < 3) {
    stop("Usage: pgs_plot.R <cnplot.txt> <output.png> <sample_id>")
}

file_in <- args[1]
file_out <- args[2]
sampleid <- args[3]

dat_nor <- read.table(file_in, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

dat_nor$chrom <- as.character(dat_nor$chrom)
dat_nor$chrom[dat_nor$chrom == "24"] <- "Y"
dat_nor$chrom[dat_nor$chrom == "23"] <- "X"
dat_nor$chrom <- sub("^chr", "", dat_nor$chrom)
dat_nor$chrom <- paste0("chr", dat_nor$chrom)
dat_nor$plot_log[dat_nor$plot_log < 0] <- 0

snp.data <- toGRanges(dat_nor[, c("chrom", "chrompos", "endpos", "plot_log")])

png(filename = file_out, width = 3000, height = 500, units = "px")
on.exit(dev.off(), add = TRUE)

kp <- plotKaryotype(
    genome = "hg19",
    plot.type = 4,
    ideogram.plotter = NULL,
    labels.plotter = NULL,
    main = sampleid,
    cex = 3
)
kpAddCytobandsAsLine(kp)
kpAddChromosomeNames(
    kp,
    chr.names <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12,
                   13, 14, 15, 16, 17, 18, 19, 20, 21, 22, "X", "Y"),
    cex = 2
)
kpPoints(
    kp,
    data = snp.data,
    y = snp.data$plot_log,
    cex = 1,
    r0 = 0,
    r1 = 1,
    ymax = 6,
    ymin = 0,
    col = colByChr(snp.data, colors = c("#458B00", "#8866DD", "#8B008B", "#FFB90F"))
)
kpAbline(kp, h = c(1, 2, 3, 4, 5, 6), r0 = 0, r1 = 1, ymax = 6, ymin = 0, col = "black")
kpAxis(kp, ymin = 0, ymax = 6, numticks = 4, cex = 2)
kpAddLabels(kp, labels = "Copy Number", srt = 90, pos = 1, label.margin = 0.03, cex = 2)
