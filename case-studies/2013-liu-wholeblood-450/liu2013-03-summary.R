library(dplyr)
library(ggplot2)
library(cowplot)

#########################################################
### Target Data: 686 450K DNA methylation samples 
###              (without known cell type proportions)
##########################################################

workingDir_liu2013 <- 
  "/users/shicks1/projects/methylCCPaper/case-studies/2013-liu-wholeblood-450"
dataPath <- "/users/shicks1/data/GEO/GSE42861"

# see liu2013-01-create-data-object.R for how 
#   * RGset_Liu2013.rds
#   * Mset_Liu2013.rds 
# were created.

counts450KHouseman <- read.csv(file.path(workingDir_liu2013, 
                          "ests/450k_minfi-estimateCellCounts.csv"), 
                          row.names = 1)

counts450KHicks <- read.csv(file.path(workingDir_liu2013, 
                          "ests/450k_methylCC-estimatecc.csv"),
                          row.names = 1)
rownames(counts450KHicks) <- rownames(counts450KHouseman)

df.Houseman <- tidyr::gather(cbind("samples" = rownames(counts450KHouseman), 
                                   data.frame(counts450KHouseman)), celltype, est, -samples)
df.Hicks <- tidyr::gather(cbind("samples" = rownames(counts450KHicks), 
                                data.frame(counts450KHicks)), celltype, est, -samples)


dfcombined <- dplyr::full_join( df.Houseman, df.Hicks, by = c("samples", "celltype"))
colnames(dfcombined) <- c("samples", "celltype", "Houseman", "methylCC")

dat_text <- as_tibble(dfcombined) %>% 
  group_by(celltype) %>% 
  mutate(sq_diff = (`Houseman` - `methylCC`)^2) %>% 
  summarize(rmse=mean(sq_diff)) %>% 
  mutate(label = paste0("RMSE: ", round(rmse,4))) %>% 
  mutate(x = rep(0, 6)) %>% 
  mutate(y = rep(0.75,6))
dat_text

gmat <- 
  ggplot(dfcombined, aes(x=Houseman, y = methylCC)) +  # 
  geom_abline(intercept = 0, slope = 1) + xlim(0,1) + ylim(0,1) +
  geom_point(size = 3, aes(color=celltype)) + facet_wrap(~celltype) + 
  theme(legend.position="top", legend.justification= "center") + 
  xlab("Houseman method") + ylab("methylCC") + 
  labs(title = "Model-based cell composition estimates from\nwhole blood samples (N=689, Li et al. 2013)", 
       color = "Cell type") + 
  scale_color_discrete(name = "Cell type", 
                       labels = c("B-cells", "CD4 T-cells", "CD8 T-cells", 
                                  "Granulocytes","Monocytes", "Natural killer cells"))

gmat <- gmat + geom_text(
  data    = dat_text,
  mapping = aes(x = x, y = y, label = label),
  hjust   = -0.1,
  vjust   = -1,
  size = 5, 
  fontface=2
)

pdf(file.path(workingDir_liu2013, 
              "figs/fig-liu2013-ests.pdf"), width = 8, height = 6)
print(gmat)
dev.off()


