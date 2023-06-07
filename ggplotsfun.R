# Before Seurat class, we need to know 2 very important R package:


#How to generate memes on R

# For windows users, you need to register your font before using it in R graphics (see discussion https://github.com/GuangchuangYu/meme/issues/1).

if (.Platform$OS.type == "windows") {
  windowsFonts(
    Impact = windowsFont("Impact"),
    Courier = windowsFont("Courier")
  )
}

# Hands on!


library(meme)
u <- system.file("angry8.jpg", package="meme")
meme(u, "code", "all the things!")

#More details on https://github.com/GuangchuangYu

###Operation phase 2 #########

# INSTALL GGCATS PACKAGE 
install.packages("remotes")
remotes::install_github("R-CoderDotCom/ggcats@main")

# LOAD LIBRARIES
library(ggplot2)
library(ggcats)
library(Ecdat)
library(tidyverse)
library(gganimate)
# VISUALIZE AVAILABLE CATS
theGrid <- expand.grid(1:5, 3:1)
dataFrame <- data.frame(x = theGrid[, 1],
                        y = theGrid[, 2],
                        image = c("shironeko", "venus", "hipster", "lil_bub", "maru", "colonel", "grumpy", "mouth", "pop", "pop_close", "pusheen", "pusheen_pc", "toast", "bongo", "nyancat"))
ggplot(dataFrame) +
  geom_cat(aes(x, y, cat = image), size = 4) +
  xlim(c(0.25, 5.5)) + 
  ylim(c(0.25, 3.5))


# TESTING PUSHEEN CAT AND MARU CAT
ggplot(mtcars) +
  geom_cat(aes(mpg, wt), cat = "pusheen", size = 4)

ggplot(mtcars) +
  geom_cat(aes(mpg, wt, size = cyl), cat = "maru")

######################################
#                                    #
#   Going further - annimated plots  #
#                                    #
######################################

# EXAMPLE GGCAT WITH GGANIMATE
library(Ecdat)
data(incomeInequality)
newData <-
  incomeInequality %>%
  select(Year, P99, median) %>%
  rename(income_median = median,
         income_99percent = P99) %>%
  pivot_longer(cols = starts_with("income"),
               names_to = "income",
               names_prefix = "income_")
newData$cat <- rep(NA, 132)
newData$cat[which(newData$income == "median")] <- rep(c("pop_close", "pop"), 33)
newData$cat[which(newData$income == "99percent")] <- "bongo"
ggplot(newData, aes(x = Year, y = value, group = income, color = income)) +
  geom_line(size = 2.5) +
  labs(title="Testing ggcats with gganimate", 
       subtitle = "Pop cat and Bongo cat") +
  geom_cat(aes(cat = cat), size = 6) +
  transition_reveal(Year)


## More resources at https://www.r-bloggers.com/2022/08/getting-fun-with-ggdogs-ggcats-and-gganimate/









