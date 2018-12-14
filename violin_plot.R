#violin Plot
setwd("D:/DATA/UC Irvine/JJEMERSON_CLASS_2018")
#res <- read.table("nolog_qn_tzb_tpm_over3000.txt",row.names=1, header = T)
library(ggplot2)
group <- read.table("violin_input.txt",header = T)
# Basic violin plot
a=group$condition
#b = log2(group$CRH+1)
b=log2(group$SOX2 +0.5)
df <- data.frame(a, b)
p1 <- ggplot(df, aes(x=a, y=b), fill=a) + 
  geom_violin(aes(fill = factor(a))) +
  #ggtitle("Dlx5") + 
  labs(title="SOX2",x="Group", y = "Log2(data)")+
  theme(text = element_text(size = 5)) + 
  geom_boxplot(aes(fill = factor(group$condition)), width = 0.05) + 
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.1, binwidth = 0.01) +
  scale_fill_brewer(palette = "Purples") +
  theme_classic()
#  geom_violin()
plot(p1)
