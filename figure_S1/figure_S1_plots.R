
library(tidyverse)
library(ggpubr)

# Host (human) lariat reads derived from lariat mapping of uninfected and IAV-infected samples of 
# wildtype and DBR1 mutant cell lines
lariat_table <- read_tsv('dbr1_mutant_patient_flu_samples_lariat_reads.txt.gz')

lariat_level <- lariat_table %>%
  filter(!(cell_line %in% c('STAT1', 'TLR3'))) %>%
  group_by(cell_line, treatment, replicate) %>%
  summarise(lariat_count = n(), total_mapped_reads = max(total_mapped_reads)) %>%
  ungroup() %>%
  mutate(lar_per_mil = lariat_count/(total_mapped_reads/1000000))

# Associate cell line labels with DBR1 mutation status
lariat_level$cell_line[which(lariat_level$cell_line %in% c('C26', 'C5', 'C1'))] <- 'WT'
lariat_level$cell_line[which(lariat_level$cell_line == 'KR')] <- 'I120T'
lariat_level$cell_line[which(lariat_level$cell_line %in% c('CCM', 'CCB'))] <- 'Y17H'
lariat_level$cell_line <- factor(lariat_level$cell_line, levels = c('WT', 'I120T', 'Y17H'))
lariat_level$treatment <- factor(lariat_level$treatment, levels = c('NS', 'IAV_MOI_5_8h', 'IAV_MOI_5_24h'))
lariat_level_mean <- lariat_level %>%
  group_by(cell_line, treatment) %>%
  summarise(lar_per_mil_mean = mean(lar_per_mil), lar_per_mil_se = sd(lar_per_mil)/sqrt(n()))

t.test(filter(lariat_level, cell_line == 'WT' & treatment == 'NS')$lar_per_mil,
       filter(lariat_level, cell_line == 'WT' & treatment == 'IAV_MOI_5_24h')$lar_per_mil)
t.test(filter(lariat_level, cell_line == 'I120T' & treatment == 'NS')$lar_per_mil,
       filter(lariat_level, cell_line == 'I120T' & treatment == 'IAV_MOI_5_24h')$lar_per_mil)
t.test(filter(lariat_level, cell_line == 'Y17H' & treatment == 'NS')$lar_per_mil,
       filter(lariat_level, cell_line == 'Y17H' & treatment == 'IAV_MOI_5_24h')$lar_per_mil)

# Plot host lariat levels for uninfected, 8hr IAV-infected, and 24hr IAV-infected samples for each cell line
ggplot(lariat_level_mean, aes(x=treatment, y=lar_per_mil_mean)) + geom_col() +
  geom_errorbar(aes(ymin=lar_per_mil_mean-lar_per_mil_se, ymax=lar_per_mil_mean+lar_per_mil_se), width=0.3) +
  facet_grid(cols=vars(cell_line)) +
  scale_x_discrete(breaks=c('NS', 'IAV_MOI_5_8h', 'IAV_MOI_5_24h'), labels=c('NI', '8h', '24h')) +
  labs(x='IAV infection', y='Human lariat reads\nper 10e6 mapped reads') +
  theme_pubr() +
  theme(axis.title = element_text(size=7), axis.text = element_text(size=7),
        strip.text = element_text(size=7))


# Host (human) lariat reads derived from lariat mapping of latent and lytic KSHV-infected samples
# GEO accession GSE280006 (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE280006)
kshv_lariat_table <- read_tsv('GSE280006_lariat_reads_merged.txt') %>%
  filter(treatment == 'DMSO')

kshv_lariat_table_summ <- kshv_lariat_table %>%
  group_by(infection, replicate) %>%
  summarise(lariat_count = n(), total_mapped_reads = max(total_mapped_reads)) %>%
  ungroup() %>%
  mutate(lariats_per_mil = lariat_count/(total_mapped_reads/1000000))

kshv_lariat_table_mean <- kshv_lariat_table_summ %>%
  group_by(infection) %>%
  summarise(lariats_per_mil_mean = mean(lariats_per_mil), lariats_per_mil_se = sd(lariats_per_mil)/sqrt(n()))

kshv_lariat_pval <- round(t.test(filter(kshv_lariat_table_summ, infection == 'LatentKSHV')$lariats_per_mil,
                           filter(kshv_lariat_table_summ, infection == 'LyticKSHV')$lariats_per_mil)$p.value, 3)

kshv_lariat_table_mean$infection <- case_match(kshv_lariat_table_mean$infection, 'LatentKSHV' ~ 'Latent', 'LyticKSHV' ~ 'Lytic')

# Plot host lariat levels for latent and lytic KSHV-infected samples
ggplot(kshv_lariat_table_mean, aes(x=infection, y=lariats_per_mil_mean)) + geom_col() +
  geom_errorbar(aes(ymin=lariats_per_mil_mean-lariats_per_mil_se, ymax=lariats_per_mil_mean+lariats_per_mil_se), width=0.4) +
  geom_bracket(xmin='Latent', xmax='Lytic', y.position=3.5, label = paste0('p = ', kshv_lariat_pval), tip.length = 0.06, label.size = 3) +
  #facet_grid(cols=vars(treatment)) +
  labs(x='KSHV infection phase', y='Lariat reads per 10e6 mapped reads') +
  theme_pubr() +
  theme(axis.text = element_text(size=8), axis.title = element_text(size=9),
        strip.text = element_text(size=9))

