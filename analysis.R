### This script analyzes the output from script_local.R. It can also easily be run on a personal computer. ###
source("script_local.R")
library(tidyMicro) #This version of tidyMicro was downloaded directly from the CharlieCarpenter/tidyMicro GitHub repository on 19 April 2021, rather than the CRAN repository
theme_set(theme_bw())

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

### Analyses not involving probiotics ###
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
# lm_bsf <- lm(FBratio ~ BSF, data = subset(phylumProp_appended.df, Treatment == "Placebo"))
# summary(lm_bsf)
# 
# lm_agebmi <- lm(BMI ~ Age, data = subset(phylumProp_appended.df, Treatment == "Placebo"))
# summary(lm_agebmi)
# 
# #T-test treatment
# t.test(FBratio ~ Treatment, data = subset(phylumProp_appended.df, Treatment != "PreTrial"))
# t.test(BSF ~ Treatment, data = subset(phylumProp_appended.df, Treatment != "PreTrial"))