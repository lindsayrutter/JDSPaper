library(ggplot2)
library(dplyr)

df <- diamonds[sample(1:nrow(diamonds), size = 1000),]
colList = scales::hue_pal()(5)
df <- df %>% select("cut","price")
ggplot(df, aes(cut, price)) + 
  geom_boxplot(aes(fill = cut)) + 
  scale_fill_manual(values = colList)


ggplot(df, aes(cut, price)) +
  #stat_boxplot(geom ='errorbar') + 
  #geom_boxplot() +
  geom_point(aes(fill = cut)) +
  scale_fill_manual(values = colList)

ggplot(df, aes(x=cut, y=price)) +
  geom_point(position=position_jitter(width=0.3), alpha=0.1, aes(fill=cut)) +
  scale_fill_manual(values=colList) +
  geom_boxplot(fill=0, outlier.size=0)

ggplot(df, aes(x=cut, y=price)) +
  stat_boxplot(geom ='errorbar') + 
  geom_boxplot(outlier.shape=NA, aes(fill=cut), alpha = 0.3) +
  geom_point(aes(fill=cut), shape=21, position=position_jitter(width=0.3), alpha=0.5) +
  scale_fill_manual(values=colList)


