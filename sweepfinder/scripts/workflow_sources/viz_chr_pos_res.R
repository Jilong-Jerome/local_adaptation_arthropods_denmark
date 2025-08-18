library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

res_file <- args[1]
spid <- args[2]
maf <- args[3]
grid <- args[4]
option <- args[5]

my_theme <- function() {
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),    # Remove major grid lines
    panel.grid.minor = element_blank(),    # Remove minor grid lines
    axis.line = element_line(color = "black"),  # Add x and y axis lines
    legend.position = "none"
  )
}

theme_set(my_theme())

## data loading
df <- read_tsv(res_file)
chroms <- df%>%
	select(chrom)%>%
	unique()
pops <- df%>%
	select(pop)%>%
	unique()
chrom_vec <- chroms$chrom
pop_vec <- pops$pop

## chromosome for all populations

for (chr in chrom_vec){
	p <- df %>%
		filter(chrom == chr)%>%
		ggplot(aes(x= location, y = LR, col = pop))+
		geom_point(size = 0.2)+
		facet_wrap(~ pop , ncol = 1,scales = "free_y")
	ggsave(paste(spid,"_maf",maf,"_w",grid,"_o",option,"_",chr,"_all_pops.pdf",sep = ""),
	p, width = 12, height = 2*length(pop_vec), limitsize = F)
}

## population for all chromosomes

for (population in pop_vec){
	p <- df %>%
		filter(pop == population)%>%
		ggplot(aes(x= location, y = LR))+
		geom_point(size = 0.2)+
		facet_wrap(~ chrom , ncol = 1, scales = "free_y")
	ggsave(paste(spid,"_maf",maf,"_w",grid,"_o",option,"_",population,"_all_chroms.pdf",sep = ""),
	p, width = 12, height = 2*length(chrom_vec), limitsize = F)
}

## per chromosome per population
for (population in pop_vec){
    for (chr in chrom_vec){
	p <- df %>%
		filter(pop == population)%>%
		filter(chrom == chr)%>%
		ggplot(aes(x= location, y = LR))+
		geom_point(size = 0.2)
	ggsave(paste(spid,"_maf",maf,"_w",grid,"_o",option,"_",population,"_",chr,"_single.pdf",sep = ""),
	p, width = 12, height = 3, limitsize = F)
    }
}
