#EDA using TidyMicro

genus_alpha <- alpha_div(tidymicro1.df, table = "Genus", min_depth = 5000)
genus_shannonH1 <- distinct(data.frame(genus_alpha1$Participant, genus_alpha1$ShannonH, genus_alpha1$Treatment))
ggplot(data = genus_shannonH, mapping = aes(x = genus_alpha.Treatment, y = genus_alpha.ShannonH)) + 
  geom_boxplot() +
  theme_bw()

tidymicro.df %>%
  ra_bars(table = "Phylum", Treatment, top_taxa = 5)

#micro_rank_sum(tidymicro.df, table = "Phylum", grp_var = Treatment)

#Phyla for the two groups
nb_phylum1 <- nb_mods(tidymicro1.df, table = "Phylum", Treatment)
micro_rocky_mtn(nb_phylum1, Treatment)
micro_forest(nb_phylum1, Treatment)

nb_phylum2 <- nb_mods(tidymicro2.df, table = "Phylum", Treatment)
micro_rocky_mtn(nb_phylum2, Treatment)
micro_forest(nb_phylum2, Treatment)

#Classes for the two groups
nb_class1 <- nb_mods(tidymicro1.df, table = "Class", Treatment)
micro_rocky_mtn(nb_class1, Treatment)
micro_forest(nb_class1, Treatment)

nb_class2 <- nb_mods(tidymicro2.df, table = "Class", Treatment)
micro_rocky_mtn(nb_class2, Treatment)
micro_forest(nb_class2, Treatment)

#Orders for the two groups
nb_order1 <- nb_mods(tidymicro1.df, table = "Order", Treatment)
micro_rocky_mtn(nb_order1, Treatment)
micro_forest(nb_order1, Treatment)

nb_order2 <- nb_mods(tidymicro2.df, table = "Order", Treatment)
micro_rocky_mtn(nb_order2, Treatment)
micro_forest(nb_order2, Treatment)

#Analyses other than probiotics
cor_rocky_mtn(subset(tidymicro.df, TestingPeriod == "1"), table = "Family", BMI)
cor_rocky_mtn(subset(tidymicro.df, TestingPeriod == "1"), table = "Family", Age)
cor_rocky_mtn(subset(tidymicro.df, TestingPeriod == "1"), table = "Family", BSF)

source("script_local.R")

#Principal component analysis
genusProp.pca <- prcomp(genusProp.df)
phylumProp.pca <- prcomp(phylumProp.df, scale. = TRUE)
autoplot(genusProp.pca, data = genusProp_appended.df, colour = 'Treatment') + theme_bw()
autoplot(phylumProp.pca, data = phylumProp_appended.df, colour = 'Treatment') + theme_bw()

#ANOVA
anova <- aov(FBratio ~ Treatment, data = phylumProp_appended.df)
TukeyHSD(anova)

#Linear regression
lm <- lm(FBratio ~ BMI + Gender + Age + Treatment + BSF, data = phylumProp_appended.df)
summary(lm)
plot(x = phylumProp_appended.df$BMI, y = phylumProp_appended.df$FBratio)

lm_bsf <- lm(FBratio ~ BSF, data = subset(phylumProp_appended.df, Treatment == "Placebo"))
summary(lm_bsf)

lm_agebmi <- lm(BMI ~ Age, data = subset(phylumProp_appended.df, Treatment == "Placebo"))
summary(lm_agebmi)

#T-test treatment
t.test(FBratio ~ Treatment, data = subset(phylumProp_appended.df, Treatment != "PreTrial"))
t.test(BSF ~ Treatment, data = subset(phylumProp_appended.df, Treatment != "PreTrial"))