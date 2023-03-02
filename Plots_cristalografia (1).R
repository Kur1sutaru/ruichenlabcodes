library(ggplot2)

tabelinhaCristal <- read.csv("~/Documents/tabelinhaCristal.csv")
mydata<-tabelinhaCristal

####P da comparacao
#Por tamanho
ggplot(mydata,
  aes_string(
    x = 'Grade',
    y = "Gene",
    size = 'P.value'
  )
) +
  geom_point(colour="red") +
ylab(NULL) + theme(axis.text.y = element_text(
  size = 12,
  angle = 0,
  hjust = 1,
  vjust = 0,
  face = "plain"
))+ scale_size(range = c(3, 8))


#Por cor
ggplot(mydata,
       aes_string(
         x = 'Grade',
         y = "Gene",
         color = 'P.value'
       )
) +
  geom_point(size=5) +
  scale_color_continuous(
    low = "red",
    high = "blue",
    name = 'P-value',
    guide = guide_colorbar(reverse = TRUE)
  ) +
  ylab(NULL) + theme(axis.text.y = element_text(
    size = 12,
    angle = 0,
    hjust = 1,
    vjust = 0,
    face = "plain"
  ))  

#########P do tecido#############
ggplot(mydata,
       aes_string(
         x = 'Grade',
         y = "Gene",
         size = 'P.value.1',
         colour="Histology"
       )
) +
  geom_point() +
  scale_colour_hue(l = 70, c = 150)+
  ylab(NULL) + theme(axis.text.y = element_text(
    size = 12,
    angle = 0,
    hjust = 1,
    vjust = 0,
    face = "plain"
  ))+ scale_size(name="P-value",range = c(3, 8))


