


library(tidyverse)
library(cowplot)
# library(dplyr)



# Normalization Example
#-------------------------------------------------------------------------------

cdata <- 
  rbind(tibble(Cell = "Cell 1", x = 1:15, y = c(.8, .3, .1, .5, .1, .3, .1, 1., .4, .4, .8, .5, .4, .2, .4)),
        tibble(Cell = "Cell 2", x = 1:15, y = c(.9, .2, .9, .4, .9, .2, .2, .9, .3, .3, .9, .6, .2, .3, .1)),
        tibble(Cell = "Cell 3", x = 1:15, y = c(.5, .4, .7, .7, .1, .5, .1, .8, .5, .4, .2, .2, .6, .4, .3)),
        tibble(Cell = "Cell 4", x = 1:15, y = c(.6, .1, .5, .3, .1, .6, .4, .9, .3, .3, .4, .7, .7, .7, .5)),
        tibble(Cell = "Cell 5", x = 1:15, y = c(.7, .2, .2, .5, .9, .7, .3, .7, .4, .2, .3, .3, .2, .4, .2)),
        tibble(Cell = "Cell 6", x = 1:15, y = c(.6, .5, .3, .2, .8, .8, .2, .8, .6, .1, .7, .7, .5, .3, .4)))

cdata$Cell = factor(cdata$Cell)
cdata %>% group_by(Cell) %>% mutate(z = as.numeric(scale(y)))


gg <- ggplot(data = cdata, aes(x=x, y= y, fill = Cell) ) +
  facet_grid(cols = vars(Cell)) +
  geom_bar(stat = "identity") + 
  theme_classic(base_size = 8) +
  xlab("Gene Expression") + ylab("Gene") +
  coord_flip() +
  theme(legend.position = "none") + 
  scale_y_continuous(n.breaks = 3)
gg
save_plot("cell_raw.png", plot = gg, base_height = 2, base_width = 7)

gg <- ggplot(data = cdata, aes(x=x, y= z, fill = Cell) ) +
  facet_grid(cols = vars(Cell)) +
  geom_bar(stat = "identity") + 
  theme_classic(base_size = 8) +
  xlab("Norm. Gene Expression") + ylab("Gene") +
  coord_flip() +
  theme(legend.position = "none") + 
  scale_y_continuous(n.breaks = 3)
gg
save_plot("cell_z.png", plot = gg, base_height = 2, base_width = 7)


# COnvert to 