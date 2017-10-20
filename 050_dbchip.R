library('DBChIP')

data("PHA4")
conds <- factor(c("emb","emb","L1", "L1"), levels=c("emb", "L1"))
path <- system.file("ext", package="DBChIP")
binding.site.list <- list()
binding.site.list[["emb"]] <- read.table(paste(path, "/emb.binding.txt", sep=""),
                                            header=TRUE)
head(binding.site.list[["emb"]])

binding.site.list[["L1"]] <- read.table(paste(path, "/L1.binding.txt", sep=""),
                                         header=TRUE)

bs.list <- read.binding.site.list(binding.site.list)

consensus.site <- site.merge(bs.list)

names(chip.data.list)                             

head(chip.data.list[["emb_rep1"]])

names(input.data.list)

dat <- load.data(chip.data.list=chip.data.list, conds=conds, consensus.site=
                        consensus.site, input.data.list=input.data.list, data.type="MCS")

dat <- get.site.count(dat)

dat <- test.diff.binding(dat)

rept <- report.peak(dat)

rept


###CAn also use BAM

library(ShortRead)
aln <- readAligned("./", pattern="emb.bam", type="Bowtie")
chip.data.list[["emb"]] <- aln
chip.data.list <- list()
chip.data.list[["emb"]] <- "/path/emb.bed.file"
chip.data.list[["L1"]] <- "/path/L1.bed.file