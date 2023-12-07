### Tutorial from https://thomaselove.github.io/431notes-2017/dataviz.html


set.seed(431001) 
library(dplyr)
library(ggplot2)
# use set.seed to ensure that we all get the same random sample 
# of 1,000 NHANES subjects in our nh_data collection

metadatainfo <- sample_n(metadatainfo, size = 63) %>%
  select(donor_ID, Age, Gender, Postmortem.time..hr., Ethnicity)

metadatainfo

ggplot(data = metadatainfo, aes(x = Age, y = Gender)) +
  geom_point()

onhinfo <- metadatainfo %>%
  filter(complete.cases(Age, Gender))

summary(onhinfo)

ggplot(data = onhinfo, aes(x = Age, y = Ethnicity, color = Gender)) +
  geom_point() +
  labs(title = "Ethnicity-Age Relationship in optic nerve samples", 
       y = "Ethnicity")


ggplot(data = onhinfo, aes(x = Age, y = Ethnicity, color = Gender)) +
  geom_point() + 
  labs(title = "Ethnicity-Age Relationship in optic nerve samples", 
       y = "Ethnicity") +
  facet_wrap(~ Gender)


ggplot(data = onhinfo, aes(x = Age, y = Ethnicity, color = Gender)) +
  geom_point() + 
  geom_smooth(method = "loess") +
  labs(title = "Ethnicity-Age Relationship in optic nerve samples", 
       y = "Ethnicity") +
  theme_bw() +
  facet_wrap(~ Gender)


ggplot(data = onhinfo, aes(x = Age)) + 
  geom_histogram() 

ggplot(data = onhinfo, aes(x = Ethnicity, y = Gender, fill = Ethnicity)) + 
  geom_boxplot() +
  labs(title = "metadata optic nerve",
       x = "Ethnicity") + 
  guides(fill = FALSE) 







