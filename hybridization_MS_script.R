library(adegenet)
library(dartR)
library(LEA)
library(vcfR)
library(poppr)

# Fri Apr 18 09:02:59 2025 ------------------------------
setwd("C:/Users/jpark107/genetics/phylogeny2/phylogeny_filter/phylogeny_all/")


setwd("C:/Users/jpark107/OneDrive - University of Tennessee/Desktop/Genetics_swo/hybridization paper")

vc<-read.vcfR("populations.snps.vcf")
g<-vcfR2genlight(vc)
g<-gl.compliance.check(g)
names<-g@ind.names
pops<-read.csv("Qbicolor_161_pops.csv", header=TRUE)
pop<-pops$pop %>% as.factor()
g@pop<-pop

#PCA plotting
library(ggplot2)
library(ggrepel)
library(dplyr)
library(cowplot) 
# Perform PCA
pca <- glPca(g, center=TRUE, scale=TRUE)
#for all samples, keep 7
barplot(pca$eig, main="eigenvalues", col=heat.colors(length(pca$eig)))
# Create a dataframe for plotting
pca_df <- data.frame(
  PC1 = pca$scores[,1],
  PC2 = pca$scores[,2],
  Population = as.factor(pop(g))  # Population as categorical variable
)
levels(pca_df[,3]) <- c(levels(pca_df[,3]), "Sample (outlier)")
# Now you can assign it
pca_df[1,3] <- "Sample (outlier)"

# Define the color palette
library('viridis')
color_palette <- viridis(29, option='turbo') 
# Select one representative per population for labeling
pca_labels <- pca_df %>%
  group_by(Population) %>%
  slice(1)  # Take the first occurrence per population

# Merge to get label positions for each point
pca_df <- pca_df %>%
  left_join(pca_labels, by = "Population", suffix = c("", "_label"))

# Variance explained (percentages for axis labels)
var_exp1 <- round(100 * pca$eig[1] / sum(pca$eig), 2)
var_exp2 <- round(100 * pca$eig[2] / sum(pca$eig), 2)
pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Population)) +
  #geom_segment(aes(x = PC1_label, y = PC2_label, xend = PC1, yend = PC2), 
  #             color = "gray50", alpha = 0.5, linetype = "dotted") +  
  geom_point(size = 4, alpha = 0.8) +  
  geom_text_repel(
    data = pca_labels,
    aes(label = Population),
    size = 5,
    family = "mono",
    max.overlaps = Inf,
    nudge_x = 0.01 * (max(pca_df$PC1) - min(pca_df$PC1)),  
    nudge_y = 0.01 * (max(pca_df$PC2) - min(pca_df$PC2)),  
    force = 2,  
    segment.size = 0.7,  
    box.padding = 1.2,
    show.legend = FALSE
  ) +
  scale_color_manual(values = color_palette) +  
  theme_minimal() +  
  theme(panel.grid.major = element_line(color = "gray85", size = 0.1),  # Lighter major grid lines
        panel.grid.minor = element_line(color = "gray92", size = 0.1))+  
  theme(
    text = element_text(family = "mono", size = 18, color='black'),  
    legend.position = "right"
  ) +
  labs(
    title = "PCA of SNP data, Dataset 3",
    x = paste0("PC1 (", var_exp1, "% variance)"),
    y = paste0("PC2 (", var_exp2, "% variance)")
  ) +
  theme(legend.title = element_blank())  

# Scree Plot (Eigenvalues Barplot)
scree_df <- data.frame(PC = 1:length(pca$eig), Eigenvalue = pca$eig)
scree_plot <- ggplot(scree_df, aes(x = PC, y = Eigenvalue)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_minimal() +
  theme(panel.grid.major = element_line(color = "gray85", size = 0),  # Lighter major grid lines
        panel.grid.minor = element_line(color = "gray92", size = 0))+  
  xlim(0,50)+
  labs(title = "Scree Plot", x = "Principal Component", y = "Eigenvalue") +
  theme(
    plot.title = element_text(size = 10, hjust = 0.5, family = "mono"),
    text = element_text(family = "mono", size = 10, color="black")
  )

# Combine PCA plot with Scree Plot in bottom corner
final_plot <- ggdraw() +
  draw_plot(pca_plot) +
  draw_plot(scree_plot, x = .1, y = 0.1, width = 0.17, height = 0.25)
final_plot


#run sNMF in LEA
setwd("C:/Users/jpark107/OneDrive - University of Tennessee/Desktop/Genetics_swo/hybridization paper/test")
gl2geno(g, outfile= "g_geno",outpath= getwd())
names<-g@ind.names
pop<-pops$pop
#detect admixture
project = snmf("g_geno.geno", K = 1:14, entropy = TRUE, repetitions = 100, seed=42, project = "new", iterations=1000000, alpha=10)
project = load.snmfProject("g_geno.snmfProject")
plot(project, col = "blue", pch = 19, cex = 1.2)
# select the best run for K = 4 clusters 
best = which.min(cross.entropy(project, K = 5)) 
my.colors <- c("tomato", "olivedrab", "gold","lightblue","pink","purple","orange") 
windows()
barchart(project, K = 5, run = best, border = NA,sort.by.Q=FALSE, space = 0, col = my.colors, xlab = "Individuals", ylab = "Ancestry proportions", main = "Ancestry matrix")-> bp 
axis(1, at = 1:length(names), labels = names, las=2, cex.axis = .5)
q_scores<-Q(project, K = 5, which.min(cross.entropy(project, K = 5)))
q<-cbind(names, pop,q_scores) %>% as.data.frame()
q$V4<-as.numeric(q$V4)
#apply threshold
q$pure<-q$V4>0.84
write.csv(q,"qscores_k5.csv")

#make txt file for STACKS populations with pure individuals
pure<-q[q$pure==TRUE,]
#remove reference samples
pure<-pure[pure$pop!='Ref. Sample-Q. bicolor',1:2]
write.table(pure, "pure.txt", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
#then run STACKS externally
stacks-2.68/populations --in-path ./phylogeny2/ --popmap ./data/popmap4_swopure111.txt --out-path ./phylogeny2/phylogeny_filter/swopure/111pops -R 0.95 --min-gt-depth 5 --min-mac 3 --genepop --structure --ordered-export --write-single-snp


#assess heterozygosity, etc. relative to Q. bicolor ancestry proportion (V4)
#get STACKS output for just the Qbicolor samples (exclude reference samples and those with higher missing data)
setwd("C:/Users/jpark107/OneDrive - University of Tennessee/Desktop/Genetics_swo/hybridization paper/test")
setwd("C:/Users/jpark107/genetics/phylogeny2/phylogeny_filter/phylogeny_all/135")
vc<-read.vcfR("populations.snps.vcf")
g<-vcfR2genlight(vc)
g<-gl.compliance.check(g)
#filter to only include SNPs with callrate of 1 (can also be done within the STACKS populations module)
g<-gl.filter.callrate(g, threshold=1)
g@pop<-pops$pop[1:135] %>% as.factor()
g@strata<-pops$pop2[1:135] %>% as.data.frame
q<-read.csv("qscores_k5.csv")

div<-gl.report.diversity(g)
heto<-gl.report.heterozygosity(g)
# Convert div$zero_D_alpha (named vector) into a data frame
div_df <- data.frame(pop = names(div$zero_D_alpha), D_alpha = as.numeric(div$zero_D_alpha))
# Merge the two tables based on the "pop" column
merged_data <- merge(div_df, heto, by = "pop")

#get Fst values and pairwise FST means
fst<-gl.fst.pop(g)
fst<-fst %>% as.matrix()
fst_no_diag <- fst
diag(fst_no_diag) <- NA
# Create an empty vector to store the means for each population
pairwise_means <- vector("numeric", length = nrow(fst))
# Loop through each population and calculate the pairwise average Fst
for (i in 1:nrow(fst)) {
  # Get the row (excluding the diagonal value for this population)
  row_values <- fst_no_diag[i, ]
  
  # Get the column (excluding the diagonal value for this population)
  col_values <- fst_no_diag[, i]
  
  # Combine both row and column values (without double-counting)
  combined_values <- c(row_values, col_values)
  
  # Remove NA values and calculate the mean
  combined_values <- combined_values[!is.na(combined_values)]
  pairwise_means[i] <- mean(combined_values)
}
# Add the results back to the matrix for easy reference
names(pairwise_means) <- rownames(fst)
pairwise_means<-as.data.frame(pairwise_means, row.names = NULL)
pairwise_means$pop<-rownames(pairwise_means)
pairwise_means
#combine with het and div table
merged_data <- merge(merged_data, pairwise_means, by = "pop")
#add a column delineating core and edge pops
merged_data$pop3<-substr(merged_data$pop, start = 1, stop = 4)
#add q scores averaged for each pop
q_135<-q[1:135,]
qavg <- q_135 %>%
  group_by(pop) %>%
  summarise(across(V4, mean, na.rm = TRUE)) %>% as.data.frame() # Compute mean for numeric columns
merged_data <- merge(merged_data, qavg, by = "pop")


#regression of FIS vs ancestry proportion
# filter out pops with less than 3 individuals
merged_data_fil<-merged_data[merged_data$n.Ind>2,]
#transform FIS to be between 0 and 1 
merged_data_fil$FIStransf<-((merged_data_fil$FIS+1)/2)
library(betareg)
mod1 <- betareg::betareg(FIStransf ~ V4, data = merged_data_fil, link = "logit")
# Extract design matrix (n × p)
s<-summary(mod1)
s$coefficients$mean[2,4]
# Extract R-squared and p-value
r2_value <- summary(mod1)$pseudo.r.squared
p_value <- summary(mod1)$coefficients$mean[2,4]  # P-value for V4
# Create text label for plot
lm_text <- paste0("Pseudo R² = ", round(r2_value, 3), 
                  "\nP-value = ", format.pval(p_value, digits = 1, eps = 0.001))
# Create scatter plot with regression line and stats annotation
new_data <- data.frame(V4 = seq(min(merged_data_fil$V4),max(merged_data_fil$V4), length.out = 100))
# Predict values and back-transform them to original scale
new_data$predicted_transf <- predict(mod1, newdata = new_data, type = "response")
new_data$predicted_original <- 2 * new_data$predicted_transf - 1  # Back-transformation
par(family="mono")
plot<-ggplot(merged_data_fil, aes(x = V4, y = FIS)) +
  geom_point(aes(color = pop)) +  # Color points by "pop"
  geom_line(data=new_data, color='black',show.legend=FALSE,size=0.5,aes(x=V4, y = predicted_original, linetype = "logit")) +
  xlab("sNMF Ancestry Proportion") +
  ylab(expression("F"[IS])) +
  theme_minimal() +
  labs(title = "Inbreeding Coefficient as Function of Ancestry Proportion", color = "Population") +
  theme(legend.position = "right", text = element_text(family = "mono", color="black")) +
  theme(panel.grid.major = element_line(color = "gray85", size = 0.05),  # Lighter major grid lines
        panel.grid.minor = element_line(color = "gray92", size = 0.05))+  
  annotate("text", 
           x = max(merged_data_fil$V4) * 0.5, 
           y = -0.1, 
           label = lm_text, 
           size = 3, family = "mono", color="black")
)
plot



#regression of average pairwise Fst values
mod <- betareg::betareg(pairwise_means ~ V4, data = merged_data_fil, link = "logit")
# Extract design matrix (n × p)
s<-summary(mod)
s$coefficients$mean[2,4]
# Extract R-squared and p-value
r2_value <- summary(mod)$pseudo.r.squared
p_value <- summary(mod)$coefficients$mean[2,4]  # P-value for V4
# Create text label for plot
lm_text <- paste0("Pseudo R² = ", round(r2_value, 3), 
                  "\nP-value = ", format.pval(p_value, digits = 3, eps = 0.001))
# Create scatter plot with regression line and stats annotation
plot<-ggplot(merged_data_fil, aes(x = V4, y = pairwise_means)) +
  geom_point(aes(color = pop)) +  # Color points by "pop"
  geom_line(color='black',show.legend=FALSE,size=0.5,aes(y = predict(mod, merged_data_fil), linetype = "logit")) +
  xlab("sNMF Ancestry Proportion") +
  ylab(expression("F"[ST])) +
  theme_minimal() +
  labs(title = "Fixation Index as Function of Ancestry Proportion", color = "Population") +
  theme(legend.position = "right", text = element_text(family = "mono", color="black")) +
  theme(panel.grid.major = element_line(color = "gray85", size = 0.05),  # Lighter major grid lines
        panel.grid.minor = element_line(color = "gray92", size = 0.05))+  
  annotate("text", 
           x = max(merged_data_fil$V4) * 0.4, 
           y = max(merged_data_fil$pairwise_means) * 0.6, 
           label = lm_text, 
           size = 3, family = "mono", color="black")
)
plot




#also get by individual observed heterozygosity
het<-gl.report.heterozygosity(g, method='ind')
colnames(het)[colnames(het) == "ind.name"] <- "names"
merged_data_ind <- merge(het, q_135, by = "names")
merged_data_ind$HetSitesHo<-(merged_data_ind$Ho)*7955
merged_data_ind$HomSitesHo<-(1-merged_data_ind$Ho)*7955
#binomial GLM model
modmod<-glm(cbind(HetSitesHo,HomSitesHo)~V4, data=merged_data_ind, family=binomial(link='logit'))
summary(modmod)
predictedmodel<-predict.glm(modmod,merged_data_ind,se.fit = T)
ci_lwr <- with(predictedmodel, plogis(fit - 1.96*se.fit))
ci_upr <- with(predictedmodel, plogis(fit + 1.96*se.fit))
merged_data_ind<-cbind(merged_data_ind,predictedmodel)
merged_data_ind$fit2<-plogis(merged_data_ind$fit)
merged_data_ind$lwr<-ci_lwr
merged_data_ind$upr<-ci_upr
#pseudo R2
r2_value<-1-modmod$deviance/modmod$null.deviance
#get pvalue
library(stats)
anova(modmod, test="Chisq")
lm_text <- paste0("\nPseudo R² = ", round(r2_value, 3),
                  "\np-value = <0.0001")
plot<-ggplot(merged_data_ind, aes(x = V4, y = Ho)) +
  geom_point(aes(color = pop)) +  # Color points by "pop"
  geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "grey", alpha = 0.15) +  # Neutral-colored error ribbon
  geom_line(aes(y = fit2), color = "black", size = .5) +  # Regression line
  xlab("sNMF Ancestry Proportion") +
  ylab(expression('H'[O])) +
  theme_minimal() +
  labs(title = "Observed Heterozygosity as Function of Ancestry Proportion", color = "Population") +
  theme(legend.position = "right", text = element_text(family = "mono", color="black")) +
  theme(panel.grid.major = element_line(color = "gray85", size = 0.05),  # Lighter major grid lines
        panel.grid.minor = element_line(color = "gray92", size = 0.05))+  
  annotate("text", 
           x = max(merged_data_ind$V4) * 0.5, 
           y = max(merged_data_ind$Ho) * 0.9, 
           label = lm_text, 
           hjust = 0, size = 3, family = "mono", color="black")
)
plot


#Dsuite
# get file from Dsuite Dinvestigate
"C:\Users\jpark107\genetics\data\dinvestigate1_output.txt"
setwd("C:/Users/jpark107/genetics/data")
# Read in lines from output file
lines <- readLines("dinvestigate1_output_2.txt")
head(lines, 20)
lines<-lines[18:length(lines)]
# Initialize list to hold parsed results
results <- list()
# We'll iterate in chunks of 7 lines (6 lines of data + 1 blank)
i <- 1
while (i <= length(lines) - 6) {
  # Skip if the current line is empty
  if (lines[i] == "") {
    i <- i + 1
    next
  }
  
  # Extract trio
  trio <- unlist(strsplit(lines[i], "\t"))
  P1 <- trio[1]
  P2 <- trio[2]
  P3 <- trio[3]
  
  # Extract values
  D       <- as.numeric(sub("D=", "", lines[i + 1]))
  f_d     <- as.numeric(sub("f_d=", "", strsplit(lines[i + 2], "\t")[[1]][1]))
  f_dM    <- as.numeric(sub("f_dM=", "", strsplit(lines[i + 3], "\t")[[1]][1]))
  ABBA_p  <- as.numeric(sub("ABBA_KSpval = ", "", lines[i + 4]))
  BABA_p  <- as.numeric(sub("BABA_KSpval = ", "", lines[i + 5]))
  
  # Save to list
  results[[length(results) + 1]] <- data.frame(
    P1 = P1, P2 = P2, P3 = P3,
    D = D, f_d = f_d, f_dM = f_dM,
    ABBA_KSpval = ABBA_p, BABA_KSpval = BABA_p,
    stringsAsFactors = FALSE
  )
  
  # Move to next block
  i <- i + 7
}
df <- do.call(rbind, results)
colnames(df)[colnames(df) == "P2"] <- "names"
colnames(pops)[colnames(pops) == "ID"] <- "names"
dfmer<-merge(df,q,by="names")
dfmer<-merge(dfmer, pops, by='names')
tail(dfmer)
plot(x=dfmer$D, y=dfmer$LAT)
plot<-ggplot(dfmer, aes(x = D, y = LAT)) +
  geom_point(aes(color = pop.y)) +  # Color points by "pop"
  xlab("sNMF Ancestry Proportion") +
  ylab(expression('D')) +
  theme_minimal() +
  labs(title = "D Stat as Function of Ancestry Proportion", color = "Population") +
  theme(legend.position = "right", text = element_text(family = "mono", color="black")) +
  theme(panel.grid.major = element_line(color = "gray85", size = 0.05),  # Lighter major grid lines
        panel.grid.minor = element_line(color = "gray92", size = 0.05))
)
plot
library(dplyr)
popD<-dfmer %>%
  group_by(pop.y) %>%
  summarize(mean_D = mean(D, na.rm = TRUE)) %>% as.data.frame()
colnames(popD)[colnames(popD) == "pop.y"] <- "pop"
merged_data<-merge(merged_data,popD, by='pop')

library(ggplot2)
library(dplyr)
# Calculate mean and SD
summary_stats <- dfmer %>%
  group_by(pop.y) %>%
  summarize(
    mean_D = mean(D, na.rm = TRUE),
    sd_D = sd(D, na.rm = TRUE)
  ) %>%
  mutate(group = ifelse(row_number() <= 9, "CORE", "EDGE"))

# Plot with color by group
plot<-ggplot(summary_stats, aes(x = pop.y, y = mean_D, fill = group)) +
  geom_col(width = 0.7) +
  geom_errorbar(aes(ymin = mean_D - sd_D, ymax = mean_D + sd_D),
                width = 0.2, color = "black") +
  scale_fill_manual(values = c("CORE" = "tomato", "EDGE" = "steelblue")) +
  labs(title = "Mean D by Population",
       x = "Population", y = "Mean D", fill = "") +
  theme(text = element_text(family = "mono")) +
  theme_minimal() +
  theme_minimal() +
  theme(
    legend.position = "right",
    text = element_text(family = "mono", color = "black"),
    panel.grid.major = element_line(color = "gray85", size = 0.05),
    panel.grid.minor = element_line(color = "gray92", size = 0.05),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
setwd("C:/Users/jpark107/OneDrive - University of Tennessee/Desktop/Genetics_swo/")
ggsave("Dstats_bypop_barplot.png",plot, width = 8, height = 6, dpi = 300, bg='white')

#plot for population level
df$pop2<-c("EDGE","EDGE","EDGE","EDGE","EDGE","EDGE","EDGE","EDGE","EDGE","EDGE","EDGE","EDGE","CORE","CORE","CORE","CORE","CORE","CORE","CORE","CORE","CORE")
library("tidyverse")
# Reshape the data to long format for plotting D, f_d, and f_dM
df_long <- df %>%
  select(P2, D, f_d, f_dM, pop2) %>%
  pivot_longer(cols = c(D, f_d, f_dM), names_to = "stat", values_to = "value")

# Plot
plot<-ggplot(df, aes(x = P2, y = D, fill = pop2)) +
  geom_col(position = "dodge", width = 0.7) +
  scale_fill_manual(values = c("CORE" = "tomato", "EDGE" = "steelblue")) +
  labs(title = "Patterson's D Statistic by Population",
       x = "Population", y = "Value", fill = "Group") +
  theme_minimal() +
  theme(
    text = element_text(family = "mono", color = "black"),
    panel.grid.major = element_line(color = "gray85", size = 0.05),
    panel.grid.minor = element_line(color = "gray92", size = 0.05),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )
plot
setwd("C:/Users/jpark107/OneDrive - University of Tennessee/Desktop/Genetics_swo/")
ggsave("Dstats_bypop_barplot2.png",plot, width = 8, height = 6, dpi = 300, bg='white')
t.test(df$D[df$pop2=="EDGE"],df$D[df$pop2=="CORE"], alternative= "greater")
range(df$D[df$pop2=="EDGE"])
range(df$D[df$pop2=="CORE"])




#load dataset only including "pure" Q. bicolor
setwd("C:/Users/jpark107/genetics/phylogeny2/phylogeny_filter/swopure/allpops/vcf")
setwd("C:/Users/jpark107/OneDrive - University of Tennessee/Desktop/Genetics_swo/hybridization paper")
vc<-read.vcfR("populations.snps.vcf")
g<-vcfR2genlight(vc)
g<-gl.compliance.check(g)
names<-g@ind.names
pure<-q[q$pure==TRUE & q$pop!='Ref. Sample-Q. bicolor',]
pure$pop3<-substr(pure$pop, start = 1, stop = 4)
g@pop<-pure$pop %>% as.factor()
strata<-pops$pop2 %>% as.data.frame()
g@strata<-strata
lat<-pops$LAT 
lon<-pops$LON
coords<-cbind(lat,lon) %>% as.data.frame()
g@other$latlon<-coords

#run IBD on edge and core together
ibd<-gl.ibd(g, coordinates="latlon")
ibd
#separate edge and core
edge<-g[g@strata$.=='EDGE']
core<-g[g@strata$.=='CORE']
#run IBD on each separately
ibd_edge<-gl.ibd(edge, coordinates="latlon")
ibd_core<-gl.ibd(core, coordinates="latlon")
amova<-gl.amova(glv)