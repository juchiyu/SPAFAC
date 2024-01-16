## code to prepare `example2_sDiSCA` dataset goes here
load("DisCA_EnFrAuthors.rda")

data.dx.unblnc$GroupDiCA <- sub("\\.", "_", data.dx.unblnc$GroupDiCA)

## replace group names ----------
## get color ------------------
group.match <- c("~Fr.18th" = "Fr($\\sim$18th)", "Fr.19th" = "Fr(19th)", "Fr.20th+" = "Fr(20th+)", "~Eng.18th"  = "Eng($\\sim$18th)",   "~UK.18th" = "UK($\\sim$18th)", "UK.19th" = "UK(19th)", "UK.20th+" = "UK(20th+)", "~Amrc.19th" = "Amrc($\\sim$19th)", "Amrc.20th+" = "Amrc(20th+)", "US.20th+" = "US(20th+)")
group.match2 <- c("~Fr_18th" = "Fr($\\sim$18th)", "Fr_19th" = "Fr(19th)", "Fr_20th+" = "Fr(20th+)", "~Eng_18th"  = "Eng($\\sim$18th)",   "~UK_18th" = "UK($\\sim$18th)", "UK_19th" = "UK(19th)", "UK_20th+" = "UK(20th+)", "~Amrc_19th" = "Amrc($\\sim$19th)", "Amrc_20th+" = "Amrc(20th+)", "US_20th+" = "US(20th+)")

## get color ------------------
col.idx <- c("Fr($\\sim$18th)" = "navy", "Fr(19th)" = "royalblue2", "Fr(20th+)" = "steelblue2", "Eng($\\sim$18th)" = "brown", "UK($\\sim$18th)" = "darkorange3", "UK(19th)" = "darkorange2","UK(20th+)" = "darkgoldenrod1", "Amrc($\\sim$19th)" = "firebrick", "Amrc(20th+)" = "palevioletred", "US(20th+)" = "palevioletred")

col.group.blnc <- as.matrix(recode(data.dx.blncC$GroupDiCA, !!!group.match) %>% recode(!!!col.idx))
rownames(col.group.blnc) <- data.dx.blncC$Authors

col.group.unblnc <- as.matrix(recode(data.dx.unblnc$GroupDiCA, !!!group.match2) %>% recode(!!!col.idx))
rownames(col.group.unblnc) <- data.dx.unblnc$Authors

## Change symbols
colnames(FrEnAuthors.unblnc)[7:9] <- c("$^\\prime$", "``\"$_{\\ll\\gg}$", "----")

data.dx.blncC$GroupDiCA.recode <- recode(data.dx.blncC$GroupDiCA, !!!group.match)
data.dx.unblnc$GroupDiCA.recode <- recode(data.dx.unblnc$GroupDiCA, !!!group.match2)

coul <- colorRampPalette(brewer.pal(8, "PiYG"))(25) # Pink (Neg) to Green (Pos)
unblnc.dx.disj <- makeNominalData(as.matrix(data.dx.unblnc$GroupDiCA.recode))
groupCntgy.unblnc <- t(unblnc.dx.disj) %*% FrEnAuthors.unblnc
rownames(groupCntgy.unblnc) <- sub(".","", rownames(groupCntgy.unblnc))


example2_sDiSCA_counts <- FrEnAuthors.unblnc
example2_sDiSCA_groups <- data.dx.unblnc$GroupDiCA.recode

usethis::use_data(example2_sDiSCA_counts, overwrite = TRUE)
usethis::use_data(example2_sDiSCA_groups, overwrite = TRUE)

