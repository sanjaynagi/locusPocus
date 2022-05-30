library(tidyverse)
library(data.table)
library(glue)
library(ggpubr)

qPCR_227 = fread("qPCR_227.csv")
qPCR_227
head(qPCR_227)



### check the correct number of samples in each gene and sample ###
table(qPCR_227$Gene)
table(qPCR_227$Strain)
table(qPCR_227$Sample)

COEAE1F = qPCR_227 %>% filter(Gene == 'COEAE1F') %>%
  select(-Well) %>% 
  group_by_at(vars(-Cq)) %>%
  summarise(Cq=mean(Cq))

EF = qPCR_227 %>% filter(Gene == 'EF') %>% 
  select(-Well) %>% 
  group_by_at(vars(-Cq)) %>%
  summarise(Cq=mean(Cq))

S7 = qPCR_227 %>% filter(Gene == 'S7') %>% 
  select(-Well) %>% 
  group_by_at(vars(-Cq)) %>%
  summarise(Cq=mean(Cq))

EF$Cq
#average reference values for each sample
refs = (S7$Cq+EF$Cq)/2

#calculate delta CT values by taking ref average and substracting Cq of 6aa1
COEAE1F$deltaCT = refs-COEAE1F$Cq

#subset to each group
COEAE1F_BAK = COEAE1F[COEAE1F$Strain == 'Bakaridjan',]
COEAE1F_KIS = COEAE1F[COEAE1F$Strain == 'Kisumu',]

#### sd and confidence interval s
errbak = sd(COEAE1F_BAK$deltaCT) *1.96 /sqrt(8)
errkis = sd(COEAE1F_KIS$deltaCT) *1.96 / sqrt(8)


#### get f 
COEAE1F_BAK_mean = mean(COEAE1F_BAK$deltaCT)
COEAE1F_KIS_mean = mean(COEAE1F_KIS$deltaCT)

### get condifence intervals
bak_lower = 2^((COEAE1F_BAK_mean - errbak) - COEAE1F_KIS_mean)
bak_upper = 2^((COEAE1F_BAK_mean + errbak) - COEAE1F_KIS_mean)

kis_lower = 2^((COEAE1F_KIS_mean - errkis) - COEAE1F_KIS_mean)
kis_upper = 2^((COEAE1F_KIS_mean + errkis) - COEAE1F_KIS_mean)


### get delta delta Ct values (log fold change to base amplication efficiency)
BAK_KIS = COEAE1F_BAK_mean-COEAE1F_KIS_mean
KIS_BAK = COEAE1F_KIS_mean-COEAE1F_BAK_mean 

KIS_KIS = COEAE1F_KIS_mean-COEAE1F_KIS_mean


## get fold change
FCBak_Kis = 2^BAK_KIS
FCBak_Kis

2^KIS_BAK


FC = data.frame("Bakaridjan"= FCBak_Kis)

FC = gather(FC, "Genotype", "Fold Change")
FC$lower = c(bak_lower)
FC$upper = c(bak_upper)

### Test for homogenous variances ####
bartlett.test(list(COEAE1F_WT$deltaCT, 
                   COEAE1F_DUP$deltaCT, 
                   COEAE1F_HET$deltaCT))

###### test for normality #### = normal
shapiro.test(COEAE1F_WT$deltaCT)
shapiro.test(COEAE1F_HET$deltaCT)
shapiro.test(COEAE1F_DUP$deltaCT)


##### perform pairwise students t test for each comparison on deltaCTs####
t.test(COEAE1F_BAK$deltaCT, COEAE1F_KIS$deltaCT)
  t.test(COEAE1F_WT$deltaCT, COEAE1F_HET$deltaCT)
t.test(COEAE1F_DUP$deltaCT, COEAE1F_HET$deltaCT)
### ggpubr t test for plotting significance *** ### same results (not needed in end)
compare_means(deltaCT ~ Strain, 
              data=COEAE1F,  method="t.test")


#### plot ####
FC$Genotype = c('Homozygous swept haplotype', 'Heterozygote', 'Homozygous wild type')
colnames(FC)[2] = "fold_change"


pdf("COEAE1F_Fold_change.pdf")
plt = FC %>% ggplot(., aes(x=reorder(Genotype, -fold_change), y=fold_change, fill=Genotype)) +
  geom_bar(stat='identity', position='dodge', color='black') +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=0.2, alpha=0.8) + 
  geom_hline(yintercept = 1, linetype='dashed', alpha=0.3)+
  ylim(c(0,3.5)) +
  theme_light() +
  theme(legend.position = 'None',
        text = element_text(size=14))+
  labs(x="Genotype", y='Relative Fold Change') +
  scale_fill_brewer(palette = "Set2") +
  geom_signif(y_position=c(3.4, 2.6, 3.1), xmin=c(0.8, 1.8, 0.8), xmax=c(3.2, 3.2, 2.2),
              annotation=c("***", "**", "NS"), tip_length=0.02)

print(plt)
dev.off() 

ggsave("COEAE1F_Fold_change.png", plt)
