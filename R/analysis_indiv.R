source("script_local.R")

#Use this script to do analysis for individuals

### Preliminary DESeq2 code that may or may not be further developed
# counts_indiv.df <- genusCount.df %>%
#   t() %>%
#   filter(isFALSE(endswith("NA", row.names())))
# 
# ddsFullCountTable <- DESeqDataSetFromMatrix( countData = counts_indiv.df,
#                                              colData = clinical.df,
#                                              design = ~ Participant + Treatment)
# 
# output <- DESeq(ddsFullCountTable)
# result <- results(output)

### Getting changes in Shannon's H at genus and phylum level

genus_shannonH <- analyzeDiversityChange(tidymicro.df, "Genus")
phylum_shannonH <- analyzeDiversityChange(tidymicro.df, "Phylum")

boxplot(genus_shannonH[[1]],
        main = "Genus-level diversity change",
        sub = paste("p = ", round(genus_shannonH[[2]]$p.value, 4)),
        xlab = "Phase",
        ylab = "Change in Shannon's H")

boxplot(phylum_shannonH[[1]],
        main = "Phylum-level diversity change",
        sub = paste("p = ", round(phylum_shannonH[[2]]$p.value, 4)),
        xlab = "Phase",
        ylab = "Change in Shannon's H")

### F/B Ratio

fb <- getFBchange(phylumProp_appended.df)
boxplot(fb[[1]],
        main = "F/B ratio change",
        sub = paste("p = ", round(fb[[2]]$p.value, 4)),
        xlab = "Phase",
        ylab = "Change in F/B ratio")

### Getting taxon-level changes

# Performing statistical analysis on proportion matrices
phylumChanges <- phylumProp_appended.df %>%
  getTaxonChange() %>%
  evalTaxonChange()
classChanges <- classProp_appended.df %>%
  getTaxonChange() %>%
  evalTaxonChange()
orderChanges <- orderProp_appended.df %>%
  getTaxonChange() %>%
  evalTaxonChange()
familyChanges <- familyProp_appended.df %>%
  getTaxonChange() %>%
  evalTaxonChange()
genusChanges <- genusProp_appended.df %>%
  getTaxonChange() %>%
  evalTaxonChange()

# Analysis on L. plantarum
Lplantarum.mat <- genusCount.df %>%
  appendData() %>%
  filterIncompletes() %>%
  findDiffs(measure = "Firmicutes/Bacilli/Lactobacillales/Lactobacillaceae/Lactiplantibacillus",
            fold = F)

# Analysis on L. plantarum persistence
Lp_present.df <- getLPstatus()
Lp.chisq <- with(Lp_present.df, chisq.test(Treatment, Lp_present))
Lp_present_count.df <- Lp_present.df %>%
  with(table(Treatment, Lp_present)) %>%
  data.frame()
ggplot(data = Lp_present_count.df, mapping = aes(x = Treatment, y = Freq, fill = Lp_present)) +
  geom_bar(stat = "identity", position = "fill") +
  ggtitle("Detection of L. plantarum by treatment",
          subtitle = paste("p =", Lp.chisq$p.value)) +
  scale_fill_viridis(discrete = TRUE)

  # Differences between study periods
Lp_present1.df <- Lp_present.df %>%
  filter(TestingPeriod %in% c("1","2")) %>%
  mutate(Treatment = ifelse(startsWith(Treatment, "Pre"), "PreTrial", Treatment))
Lp1.chisq <- with(Lp_present1.df, chisq.test(Treatment, Lp_present))
Lp_present_count1.df <- Lp_present1.df %>%
  with(table(Treatment, Lp_present)) %>%
  data.frame()
ggplot(data = Lp_present_count1.df, mapping = aes(x = Treatment, y = Freq, fill = Lp_present)) +
  geom_bar(stat = "identity", position = "fill") +
  ggtitle("Detection of L. plantarum by treatment",
          subtitle = "Study Period 1") +
  scale_fill_viridis(discrete = TRUE)

Lp_present2.df <- Lp_present.df %>%
  filter(TestingPeriod %in% c("3","4")) %>%
  mutate(Treatment = ifelse(startsWith(Treatment, "Pre"), "PreTrial", Treatment))
Lp2.chisq <- with(Lp_present2.df, chisq.test(Treatment, Lp_present))
Lp_present_count2.df <- Lp_present2.df %>%
  with(table(Treatment, Lp_present)) %>%
  data.frame()
ggplot(data = Lp_present_count2.df, mapping = aes(x = Treatment, y = Freq, fill = Lp_present)) +
  geom_bar(stat = "identity", position = "fill") +
  ggtitle("Detection of L. plantarum by treatment",
          subtitle = "Study Period 2") +
  scale_fill_viridis(discrete = TRUE)

# Tests on other variables related to persistence
Lp_present_active.df <- filter(Lp_present.df, Treatment == "Lp299v")
with(Lp_present_active.df, chisq.test(Overweight, Lp_present))
with(Lp_present_active.df, chisq.test(Site, Lp_present))
with(Lp_present_active.df, chisq.test(Gender, Lp_present))
t.test(BSF ~ Lp_present, data = Lp_present_active.df)



# Principal component analysis
## For now, this has been superseded by RDA
# genusProp_appended_start.df <- filter(genusProp_appended.df, Treatment == "PreTrial_1")
# genusProp_appended_Lp299v.df <- filter(genusProp_appended.df, Treatment == "Lp299v")
# 
# genusProp_appended_start_pca.df <- genusProp_appended_start.df %>%
#   select(-c(colnames(clinical.df[,-1]),
#             Lp_present,
#             Firmicutes.Bacilli.Lactobacillales.Lactobacillaceae.Lactiplantibacillus)) %>%
#   select_if(function(col) max(col) != 0)
# 
# genusProp_appended_Lp299v_pca.df <- genusProp_appended_Lp299v.df %>%
#   select(-c(colnames(clinical.df[,-1]),
#             Lp_present,
#             Firmicutes.Bacilli.Lactobacillales.Lactobacillaceae.Lactiplantibacillus)) %>%
#   select_if(function(col) max(col) != 0)
# 
# genus_Lp.pca <- prcomp(genusProp_appended_Lp299v_pca.df,
#                        scale = TRUE)
# autoplot(genus_Lp.pca, data = genusProp_appended_Lp299v.df, colour = 'Lp_present') +
#   ggtitle("Genus-level scaled PCA",
#           subtitle = "Treatment phase") +
#   scale_color_viridis(discrete = TRUE)
# 
# genus_start.pca <- prcomp(genusProp_appended_start_pca.df,
#                           scale = TRUE)
# autoplot(genus_start.pca, data = genusProp_appended_start.df, colour = 'BMI') +
#   ggtitle("Genus-level PCA",
#           subtitle = "Testing Period 1") +
#   scale_color_viridis(discrete = FALSE)
# autoplot(genus_start.pca, data = genusProp_appended_start.df, colour = 'Age') +
#   ggtitle("Genus-level PCA") +
#   scale_color_viridis(discrete = FALSE)
# autoplot(genus_start.pca, data = genusProp_appended_start.df, colour = 'Gender') +
#   ggtitle("Genus-level scaled PCA",
#           subtitle = "Testing Period 1") +
#   scale_color_viridis(discrete = TRUE)
