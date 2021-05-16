### This script analyzes the output from script_local.R. It can also easily be run on a personal computer. ###
source("script_local.R")

### Diversity analysis ###

#Alpha diversity

genus_alpha1 <- alpha_div(tidymicro1.df, table = "Genus", min_depth = 5000)
genus_shannonH1 <- distinct(data.frame(genus_alpha1$Participant, genus_alpha1$ShannonH, genus_alpha1$Treatment))
ggplot(data = genus_shannonH1, mapping = aes(x = genus_alpha1.Treatment, y = genus_alpha1.ShannonH)) + 
  geom_boxplot() +
  ggtitle("Boxplot of Shannon's H in different treatments",
          subtitle = "Period 1, Genus-level analysis")

genus_alpha2 <- alpha_div(tidymicro2.df, table = "Genus", min_depth = 5000)
genus_shannonH2 <- distinct(data.frame(genus_alpha2$Participant, genus_alpha2$ShannonH, genus_alpha2$Treatment))
ggplot(data = genus_shannonH2, mapping = aes(x = genus_alpha2.Treatment, y = genus_alpha2.ShannonH)) + 
  geom_boxplot() +
  ggtitle("Boxplot of Shannon's H in different treatments",
          subtitle = "Period 2, Genus-level analysis")

#Taxon barplots

tidymicro1.df %>%
  ra_bars(table = "Order", Treatment, top_taxa = 10) +
  ggtitle("Stacked bar chart of taxa in different treatments",
          subtitle = "Period 1, Order-level analysis")

tidymicro2.df %>%
  ra_bars(table = "Order", Treatment, top_taxa = 10) +
  ggtitle("Stacked bar chart of taxa in different treatments",
          subtitle = "Period 2, Order-level analysis")

### Forest plots and Rocky Mountain plots ###

#Phylum-level analysis
nb_phylum1 <- nb_mods(tidymicro1.df, table = "Phylum", Treatment)
micro_rocky_mtn(nb_phylum1, Treatment) + ggtitle("Rocky Mountain plot of Lp299v strain impact",
                                                 subtitle = "Period 1, Phylum-level analysis")
micro_forest(nb_phylum1, Treatment) + ggtitle("Forest plot of Lp299v strain impact",
                                              subtitle = "Period 1, Phylum-level analysis")

nb_phylum2 <- nb_mods(tidymicro2.df, table = "Phylum", Treatment)
micro_rocky_mtn(nb_phylum2, Treatment) + ggtitle("Rocky Mountain plot of Lp299v strain impact",
                                                 subtitle = "Period 2, Phylum-level analysis")
micro_forest(nb_phylum2, Treatment) + ggtitle("Forest plot of Lp299v strain impact",
                                              subtitle = "Period 2, Phylum-level analysis")

#Class-level analysis
nb_class1 <- nb_mods(tidymicro1.df, table = "Class", Treatment)
micro_rocky_mtn(nb_class1, Treatment) + ggtitle("Rocky Mountain plot of Lp299v strain impact",
                                                subtitle = "Period 1, Class-level analysis")
micro_forest(nb_class1, Treatment) + ggtitle("Forest plot of Lp299v strain impact",
                                             subtitle = "Period 1, Class-level analysis")

nb_class2 <- nb_mods(tidymicro2.df, table = "Class", Treatment)
micro_rocky_mtn(nb_class2, Treatment) + ggtitle("Rocky Mountain plot of Lp299v strain impact",
                                                subtitle = "Period 2, Class-level analysis")
micro_forest(nb_class2, Treatment) + ggtitle("Forest plot of Lp299v strain impact",
                                             subtitle = "Period 2, Class-level analysis")

#Order-level analysis
nb_order1 <- nb_mods(tidymicro1.df, table = "Order", Treatment)
micro_rocky_mtn(nb_order1, Treatment, alpha = 0.001) + ggtitle("Rocky Mountain plot of Lp299v strain impact",
                                                               subtitle = "Period 1, Order-level analysis")
micro_forest(nb_order1, Treatment) + ggtitle("Forest plot of Lp299v strain impact",
                                             subtitle = "Period 1, Order-level analysis")


nb_order2 <- nb_mods(tidymicro2.df, table = "Order", Treatment)
micro_rocky_mtn(nb_order2, Treatment, alpha = 0.001) + ggtitle("Rocky Mountain plot of Lp299v strain impact",
                                                               subtitle = "Period 2,Order-level analysis")
micro_forest(nb_order2, Treatment) + ggtitle("Forest plot of Lp299v strain impact",
                                             subtitle = "Period 2, Order-level analysis")

### Redundancy analysis ###

# Comparing baselines

genusProp_appended_baselines_pca.df <- filter(genusProp_appended.df, Treatment %in% c("PreTrial_2", "PreTrial_1")) %>%
  select(-c(colnames(clinical.df[,-1])))
clinical_baselines.df <- filter(clinical.df, Treatment %in% c("PreTrial_2", "PreTrial_1"))

baselines.rda <- rda(genusProp_appended_baselines_pca.df ~ Treatment, data = clinical_baselines.df)
baselines_rda.anova <- anova(baselines.rda)

# Gender

genusProp_appended_start.df <- filter(genusProp_appended.df, Treatment == "PreTrial_1")

genusProp_appended_start_pca.df <- genusProp_appended_start.df %>%
  select(-c(colnames(clinical.df[,-1]),
            Lp_present,
            Firmicutes.Bacilli.Lactobacillales.Lactobacillaceae.Lactiplantibacillus)) %>%
  select_if(function(col) max(col) != 0)

gender.rda <- rda(genusProp_appended_start_pca.df ~ Gender,
                  data = genusProp_appended_start.df,
                  scale = TRUE)

RsquareAdj(gender.rda)
anova(gender.rda)

# Lp299v
clinical_Lp299v.df <- filter(clinical.df, Treatment == "Lp299v")

genusProp_appended_Lp299v_pca.df <- genusCount_appended.df %>%
  filter(Treatment == "Lp299v") %>%
  select(-c(colnames(clinical.df[,-1]),
            Lp_present)) %>%
  select_if(function(col) max(col) != 0) %>%
  select(!ends_with("Lactiplantibacillus")) %>%
  countToProp() #converting back to prop table after eliminating L. plantarum just to control for ANY extraneous impact of including Lp in calculations

lp299v.rda <- rda(genusProp_appended_Lp299v_pca.df ~ Lp_present,
                  scale = TRUE,
                  data = clinical_Lp299v.df,
                  na.rm = TRUE)

anova <- anova(lp299v.rda, permutations = 9999)
RsquareAdj(lp299v.rda)
lp299v_rda_plot <- autoplot(lp299v.rda, geom = "point")
lp299v_rda_plot$layers[[1]]$geom$default_aes$size <- 1.5
lp299v_rda_plot$layers[[1]]$geom$default_aes$fill <- "#FFFFFF"
lp299v_rda_plot

#Note that ANOVA is calculated stochastically for RDA objects. Therefore, the one you calculate may vary slightly from the reported value.
plot(lp299v.rda)
title("Redundancy analysis on Lp299v treatment phase (Genus-level analysis)",
      sub = paste("p =", anova$`Pr(>F)`[1]))

# Lp299v/baseline
clinical_Lp_baseline.df <- filter(clinical.df, (Treatment == "PreTrial_1" & TestingOrder == "Lp299v - Placebo") |
                                    (Treatment == "PreTrial_2" & TestingOrder == "Placebo - Lp299v"))

genusProp_appended_Lp_baseline_pca.df <- genusCount_appended.df %>%
  filter((Treatment == "PreTrial_1" & TestingOrder == "Lp299v - Placebo") |
           (Treatment == "PreTrial_2" & TestingOrder == "Placebo - Lp299v")) %>%
  select(-c(colnames(clinical.df[,-1]),
            Lp_present)) %>%
  select_if(function(col) max(col) != 0) %>%
  select(!ends_with("Lactiplantibacillus")) %>%
  countToProp() #converting back to prop table after eliminating L. plantarum just to control for ANY extraneous impact of including Lp in calculations

lp_baseline.rda <- rda(genusProp_appended_Lp_baseline_pca.df ~ Lp_present_treat,
                  scale = TRUE,
                  data = clinical_Lp_baseline.df,
                  na.rm = TRUE)

anova <- anova(lp_baseline.rda, permutations = 9999)
RsquareAdj(lp_baseline.rda)

#Note that ANOVA is calculated stochastically for RDA objects. Therefore, the one you calculate may vary slightly from the reported value.
plot(lp_baseline.rda)
title("Redundancy analysis on baseline (Genus-level analysis)",
      sub = paste("p =", anova$`Pr(>F)`[1]))

### Extraneous tidyMicro analyses involving probiotics ###
tidymicro_Lp299v.df <- filter(tidymicro.df, Treatment == "Lp299v")%>%
  mutate(Overweight = factor(ifelse(Overweight, "yes", "no"))) %>%
  mutate(Lp_present = factor(ifelse(Lp_present, "yes", "no"))) %>%
  mutate(Gender = factor(Gender))

# Correlation matrix

cor_heatmap(tidymicro_Lp299v.df, table = "Order", BMI, Age, BSF) +
  scale_fill_viridis(option = "magma") +
  ggtitle("Order-level Correlation Plot",
          subtitle = "Lp299v phase")

# Beta diversity
Lp299v_beta <- beta_div(tidymicro_Lp299v.df, table = "Genus", method = "bray")
Lp299v_beta %>%
  beta_heatmap(micro_set = tidymicro_Lp299v.df, Lp_present) +
  scale_fill_viridis(option = "magma") +
  ggtitle("Beta diversity analysis by Lp299v status",
          subtitle = "Genus-level analysis")

# Taxa barplot
ra_bars(tidymicro_Lp299v.df,
        table = "Genus",
        Lp_present,
        top_taxa = 15)

# Taxa barplot (comparing baselines)
tidymicro_Lp_baseline.df <- filter(tidymicro.df, (Treatment == "PreTrial_1" & TestingOrder == "Lp299v - Placebo") |
                                     (Treatment == "PreTrial_2" & TestingOrder == "Placebo - Lp299v"))

ra_bars(tidymicro_Lp_baseline.df,
        table = "Genus",
        Lp_present_treat,
        top_taxa = 15)

### Analyses not involving probiotics ###
tidymicro_period1.df <- filter(tidymicro.df, TestingPeriod == "1")%>%
  mutate(Overweight = factor(ifelse(Overweight, "yes", "no"))) %>%
  mutate(Gender = factor(Gender))

cor_heatmap(tidymicro_period1.df, table = "Order", BMI, Age, BSF) +
  scale_fill_viridis(option = "magma") +
  ggtitle("Order-level Correlation Plot",
          subtitle = "Testing Period 1")

# Alpha diversity
period1_alpha <- alpha_div(tidymicro_period1.df,
                           table = "Genus",
                           min_depth = 5000,
                           min_goods = 80)
alphareg <- micro_alpha_reg(period1_alpha,
                            table = "Genus",
                            Age,
                            BMI,
                            BSF)

# Beta diversity
period1_beta <- beta_div(tidymicro_period1.df, table = "Genus", method = "bray")
period1_beta %>%
  beta_heatmap(micro_set = tidymicro_period1.df, Gender) +
  scale_fill_viridis(option = "magma") +
  ggtitle("Beta diversity analysis by gender",
          subtitle = "Genus-level analysis, Period 1")

cor_rocky_mtn(subset(tidymicro.df, TestingPeriod == "1"), table = "Family", BMI)
nb_class_bmi <- nb_mods(subset(tidymicro.df, TestingPeriod == "1"), table = "Class", Overweight)
micro_rocky_mtn(nb_class_bmi, Overweight, alpha = 0.001)
cor_rocky_mtn(subset(tidymicro.df, TestingPeriod == "1"), table = "Family", Age)
cor_rocky_mtn(subset(tidymicro.df, TestingPeriod == "1"), table = "Family", BSF)

### Old code that may be recycled ###

# #Principal component analysis
# genusProp.pca <- prcomp(genusProp.df)
# phylumProp.pca <- prcomp(phylumProp.df, scale. = TRUE)
# autoplot(genusProp.pca, data = genusProp_appended.df, colour = 'Treatment') + theme_bw()
# autoplot(phylumProp.pca, data = phylumProp_appended.df, colour = 'Treatment') + theme_bw()
# 
# #ANOVA
# anova <- aov(FBratio ~ Treatment, data = phylumProp_appended.df)
# TukeyHSD(anova)
# 
# #Linear regression
# lm <- lm(FBratio ~ BMI + Gender + Age + Treatment + BSF, data = phylumProp_appended.df)
# summary(lm)
# plot(x = phylumProp_appended.df$BMI, y = phylumProp_appended.df$FBratio)
# 
# lm_bsf <- lm(FBratio ~ BSF, data = subset(phylumProp_appended.df, Treatment == "PreTrial_1"))
# summary(lm_bsf)
# 
# lm_agebmi <- lm(BMI ~ Age, data = subset(phylumProp_appended.df, Treatment == "Placebo"))
# summary(lm_agebmi)
# 
# #T-test treatment
# t.test(FBratio ~ Treatment, data = subset(phylumProp_appended.df, Treatment != "PreTrial"))
# t.test(BSF ~ Treatment, data = subset(phylumProp_appended.df, Treatment != "PreTrial"))