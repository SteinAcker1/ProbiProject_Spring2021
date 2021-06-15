source("script_local.R")

### Basic demographic info (Table 1 data) ###

# Age by group

"Age"
clinical.df %>%
  filter(Treatment == "PreTrial_1") %>%
  with(tapply(Age, TestingOrder, summary))

# BMI by group

"BMI"
clinical.df %>%
  filter(Treatment == "PreTrial_1") %>%
  with(tapply(BMI, TestingOrder, summary))

# BSF by group and treatment

"BSF at Lp299v baseline"
clinical.df %>%
  filter((Treatment == "PreTrial_1" & TestingOrder == "Lp299v - Placebo") |
           (Treatment == "PreTrial_2" & TestingOrder == "Placebo - Lp299v")) %>%
  with(tapply(BSF, TestingOrder, summary)) #Pre-Lp299v

"BSF at placebo baseline"
clinical.df %>%
  filter((Treatment == "PreTrial_2" & TestingOrder == "Lp299v - Placebo") |
           (Treatment == "PreTrial_1" & TestingOrder == "Placebo - Lp299v")) %>%
  with(tapply(BSF, TestingOrder, summary)) #Pre-Placebo

"BSF after Lp299v treatment"
clinical.df %>%
  filter(Treatment == "Lp299v") %>%
  with(tapply(BSF, TestingOrder, summary)) #Lp299v

"BSF after placebo treatment"
clinical.df %>%
  filter(Treatment == "Placebo") %>%
  with(tapply(BSF, TestingOrder, summary)) #Placebo

# Diversity

genus_alpha.df <- tidymicro.df %>%
  alpha_div(table = "Genus", min_depth = 5000) %>%
  select(c("Treatment", "TestingOrder", "ShannonH", "Participant", "Lp_present_treat")) %>%
  distinct()

"Shannon's H at Lp299v baseline"
genus_alpha.df %>%
  filter((Treatment == "PreTrial_1" & TestingOrder == "Lp299v - Placebo") |
           (Treatment == "PreTrial_2" & TestingOrder == "Placebo - Lp299v")) %>%
  with(tapply(ShannonH, TestingOrder, summary)) #Pre-Lp299v

"Shannon's H at placebo baseline"
genus_alpha.df %>%
  filter((Treatment == "PreTrial_2" & TestingOrder == "Lp299v - Placebo") |
           (Treatment == "PreTrial_1" & TestingOrder == "Placebo - Lp299v")) %>%
  with(tapply(ShannonH, TestingOrder, summary)) #Pre-Placebo

"Shannon's H after Lp299v treatment"
genus_alpha.df %>%
  filter(Treatment == "Lp299v") %>%
  with(tapply(ShannonH, TestingOrder, summary)) #Lp299v

"Shannon's H after placebo treatment"
genus_alpha.df %>%
  filter(Treatment == "Placebo") %>%
  with(tapply(ShannonH, TestingOrder, summary)) #Placebo

"Lactiplantibacillus in gut during treatment - Pre-Lp299v Shannon's H"
genus_alpha.df %>%
  filter((Treatment == "PreTrial_1" & TestingOrder == "Lp299v - Placebo") |
           (Treatment == "PreTrial_2" & TestingOrder == "Placebo - Lp299v")) %>%
  with(tapply(ShannonH, Lp_present_treat, summary)) #Pre-Lp299v, based on active status

"Lactiplantibacillus in gut during treatment - Post-Lp299v Shannon's H"
genus_alpha.df %>%
  filter(Treatment == "Lp299v") %>%
  with(tapply(ShannonH, Lp_present_treat, summary)) #Lp299v, based on active status

# F/B Ratio

fb.df <- phylumProp_appended.df %>%
  mutate(FBratio = Firmicutes / Bacteroidota)

"F/B ratio at Lp299v baseline"
fb.df %>%
  filter((Treatment == "PreTrial_1" & TestingOrder == "Lp299v - Placebo") |
           (Treatment == "PreTrial_2" & TestingOrder == "Placebo - Lp299v")) %>%
  with(tapply(FBratio, TestingOrder, summary)) #Pre-Lp299v

"F/B ratio at placebo baseline"
fb.df %>%
  filter((Treatment == "PreTrial_2" & TestingOrder == "Lp299v - Placebo") |
           (Treatment == "PreTrial_1" & TestingOrder == "Placebo - Lp299v")) %>%
  with(tapply(FBratio, TestingOrder, summary)) #Pre-Placebo

"F/B ratio after Lp299v treatment"
fb.df %>%
  filter(Treatment == "Lp299v") %>%
  with(tapply(FBratio, TestingOrder, summary)) #Lp299v

"F/B ratio after placebo treatment"
fb.df %>%
  filter(Treatment == "Placebo") %>%
  with(tapply(FBratio, TestingOrder, summary)) #Placebo

### P-values: t-tests ###
"Difference between placebo and Lp299v - t-test, Shannon's H"
genus_alpha.df %>%
  filter(Treatment %in% c("Placebo", "Lp299v")) %>%
  with(t.test(ShannonH ~ Treatment)) # based on treatment

"Difference between baselines - t-test, Shannon's H"
genus_alpha.df %>%
  filter(Treatment %in% c("PreTrial_1", "PreTrial_2")) %>%
  with(t.test(ShannonH ~ Treatment)) # based on baselines

"Lactiplantibacillus in gut during treatment - Pre-Lp299v t-test, Shannon's H"
genus_alpha.df %>%
  filter((Treatment == "PreTrial_1" & TestingOrder == "Lp299v - Placebo") |
           (Treatment == "PreTrial_2" & TestingOrder == "Placebo - Lp299v")) %>%
  with(t.test(ShannonH ~ Lp_present_treat)) # based on baselines to Lp299v status

"Lactiplantibacillus in gut during treatment - Post-Lp299v t-test, Shannon's H"
genus_alpha.df %>%
  filter(Treatment == "Lp299v") %>%
  with(t.test(ShannonH ~ Lp_present_treat)) # based on Lp299v status

### P-values: RDA ###

"Lactiplantibacillus in gut during treatment - Pre-Lp299v RDA"

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

anova(lp_baseline.rda, permutations = 9999)

"Lactiplantibacillus in gut during treatment - Post-Lp299v RDA"

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

anova(lp299v.rda, permutations = 9999)

### Figures ###

# Fig. 2

Lp_present.df <- getLPstatus()

Lp_present_count.df <- Lp_present.df %>%
  with(table(Treatment, Lp_present, TestingOrder)) %>%
  data.frame()

Fig2.df <- Lp_present_count.df %>%
  group_by(Treatment, TestingOrder) %>%
  mutate(pct = 100 * Freq / sum(Freq)) %>%
  filter(Lp_present == TRUE) %>%
  select(-Lp_present)

jpeg(filename = "../output/Fig2.jpeg", width = 667, height = 453)
ggplot(data = Fig2.df, mapping = aes(fill = TestingOrder, y = pct, x = Treatment)) +
  geom_bar(stat = "identity", position = "dodge")
dev.off()

# Fig. 3A

rmplot_Lp_present <- genusProp_appended.df %>%
  getDemographicDiff(var = "Lp_present", treat = "Lp299v") %>%
  filter(str_detect(taxon, "Lactiplantibacillus$", negate = TRUE)) %>%
  shiftRockyMtn()

jpeg(filename = "../output/Fig3A.jpeg", width = 667, height = 453)
genusRockymtnLp299vpresence <- rmplot_Lp_present[[2]] +
  ylab(TeX(r"($log(1/p) * sgn(\Delta)$)")) +
  labs(tag = "A")
dev.off()

# Fig. 3B

rmplot_Lp_present_baseline <- genusProp_appended.df %>%
  mutate(Treatment = ifelse((Treatment == "PreTrial_1" & TestingOrder == "Lp299v - Placebo") |
                              (Treatment == "PreTrial_2" & TestingOrder == "Placebo - Lp299v"),
                            "PreLp299v",
                            Treatment)) %>%
  getDemographicDiff(var = "Lp_present_treat", treat = "PreLp299v") %>%
  filter(str_detect(taxon, "Lactiplantibacillus$", negate = TRUE)) %>%
  shiftRockyMtn()

jpeg(filename = "../output/Fig3B.jpeg", width = 667, height = 453)
genusRockymtnLp299vpresence_baseline <- rmplot_Lp_present_baseline[[2]] +
  ylab(TeX(r"($log(1/p) * sgn(\Delta)$)")) +
  labs(tag = "B")
dev.off()

# Fig. 4
genus_alpha_boxplot.df <- genus_alpha.df %>%
  filter((Treatment == "PreTrial_1" & TestingOrder == "Lp299v - Placebo") |
           (Treatment == "PreTrial_2" & TestingOrder == "Placebo - Lp299v") |
           Treatment == "Lp299v") %>%
  mutate(Treatment = ifelse(Treatment == "Lp299v", Treatment, "Baseline")) %>%
  filter(!is.na(Lp_present_treat))

jpeg(filename = "../output/Fig4.jpeg", width = 667, height = 453)
ggplot(data = genus_alpha_boxplot.df, mapping = aes(x = Treatment, y = ShannonH, fill = Lp_present_treat)) +
  geom_boxplot()
dev.off()

# Supp. fig. 1A

jpeg(filename = "../output/suppFig1A.jpeg", width = 667, height = 453)
ra_bars(tidymicro_Lp299v.df,
        table = "Genus",
        Lp_present_treat,
        top_taxa = 15) +
  labs(tag = "A")
dev.off()

# Supp. fig. 1B

jpeg(filename = "../output/suppFig1B.jpeg", width = 667, height = 453)
ra_bars(tidymicro_Lp_baseline.df,
        table = "Genus",
        Lp_present_treat,
        top_taxa = 15) +
  labs(tag = "B")
dev.off()

# Supp. fig. 2

jpeg(filename = "../output/suppFig2.jpeg", width = 667, height = 453)
Lp299v_beta <- beta_div(tidymicro_Lp299v.df, table = "Genus", method = "bray")
Lp299v_beta %>%
  beta_heatmap(micro_set = tidymicro_Lp299v.df, Lp_present) +
  scale_fill_viridis(option = "magma")
dev.off()

# Outputting clinical data TSV
write.table(clinical.df, file = "../output/ProGastro17_demographics.tsv", sep = "\t", row.names = FALSE)
