setwd("D:/baylor/DRG atlas/DEG")

library(dplyr)
library(readr)

#import and merge all three CSV files into one data frame
df <- list.files(path='D:/baylor/DRG atlas/DEG') %>% 
  lapply(read_csv) %>% 
  bind_rows 

# Save the merged table
write.csv(df, "degs_all_DRG.csv")
