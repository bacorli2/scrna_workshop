


library(tidyverse)
library(cowplot)
# library(dplyr)


# Set wd to base of workshop repository
here::i_am("README.md")

# Normalization Data Example
#-------------------------------------------------------------------------------
set.seed(1)
cdata <- 
  rbind(tibble(Cell = "Cell 1", x = 1:15, y = runif(15, min = 0, max = 30)),
        tibble(Cell = "Cell 2", x = 1:15, y = runif(15, min = 0, max = 30)),
        tibble(Cell = "Cell 3", x = 1:15, y = runif(15, min = 0, max = 30)),
        tibble(Cell = "Cell 4", x = 1:15, y = runif(15, min = 0, max = 30)),
        tibble(Cell = "Cell 5", x = 1:15, y = runif(15, min = 0, max = 30)),
        tibble(Cell = "Cell 6", x = 1:15, y = runif(15, min = 0, max = 30)))

cdata$Cell = factor(cdata$Cell)
cdata <- cdata %>% group_by(Cell) %>% mutate(z = log1p(10000 * y / (sum(y)*10000)))


gg <- ggplot(data = cdata, aes(x=x, y= y, fill = Cell) ) +
  facet_grid(cols = vars(Cell)) +
  geom_bar(stat = "identity") + 
  theme_classic(base_size = 8) +
  xlab("Gene ID") + ylab("Count") +
  coord_flip() +
  theme(legend.position = "none") + 
  scale_y_continuous(n.breaks = 3)
gg
save_plot(here::here("_temp_out","cell_raw_count.png"), plot = gg, base_height = 2, base_width = 7)

gg <- ggplot(data = cdata, aes(x=x, y= z, fill = Cell) ) +
  facet_grid(cols = vars(Cell)) +
  geom_bar(stat = "identity") + 
  theme_classic(base_size = 8) +
  xlab("Gene ID") + ylab("Norm. Count") +
  coord_flip() +
  theme(legend.position = "none") + 
  scale_y_continuous(n.breaks = 3)
gg
save_plot(here::here("_temp_out","cell_norm_count.png"), plot = gg, base_height = 2, base_width = 7)



# Scaling Data Example
#-------------------------------------------------------------------------------
cdata <- 
  rbind(tibble(Cell = "Cell 1", Gene = 1:15, y = c(.8, .3, .1, .5, .1, .3, .1, 1., .4, .4, .8, .5, .4, .2, .4)),
        tibble(Cell = "Cell 2", Gene = 1:15, y = c(.9, .2, .9, .4, .9, .2, .2, .9, .3, .3, .9, .6, .2, .3, .1)),
        tibble(Cell = "Cell 3", Gene = 1:15, y = c(.5, .4, .7, .7, .1, .5, .1, .8, .5, .4, .2, .2, .6, .4, .3)),
        tibble(Cell = "Cell 4", Gene = 1:15, y = c(.6, .1, .5, .3, .1, .6, .4, .9, .3, .3, .4, .7, .7, .7, .5)),
        tibble(Cell = "Cell 5", Gene = 1:15, y = c(.7, .2, .2, .5, .9, .7, .3, .7, .4, .2, .3, .3, .2, .4, .2)),
        tibble(Cell = "Cell 6", Gene = 1:15, y = c(.6, .5, .3, .2, .8, .8, .2, .8, .6, .1, .7, .7, .5, .3, .4)))

# Z_score when measurements are grouped by gene
cdata$Cell = factor(cdata$Cell)
cdata$Gene = factor(cdata$Gene)

cdata <- cdata %>% group_by(Gene) %>% mutate(z = as.numeric(scale(y)))


gg <- ggplot(data = cdata, aes(x = Gene, y= y, fill = Cell) ) +
  facet_grid(cols = vars(Cell)) +
  geom_bar(stat = "identity") + 
  theme_classic(base_size = 8) +
  ylab("Gene Expression") + xlab("Gene ID") +
  coord_flip() +
  theme(legend.position = "none") + 
  scale_y_continuous(n.breaks = 3)
gg
save_plot(here::here("_temp_out","cell_raw_scale.png"), plot = gg, base_height = 2, base_width = 7)

gg <- ggplot(data = cdata, aes(x = Gene, y = z, fill = Cell) ) +
  facet_grid(cols = vars(Cell)) +
  geom_bar(stat = "identity") + 
  theme_classic(base_size = 8) +
  ylab("Scaled. Gene Expression") + xlab("Gene ID") +
  coord_flip() +
  theme(legend.position = "none") + 
  scale_y_continuous(n.breaks = 3)
gg
save_plot(here::here("_temp_out","cell_scaled_scaled.png"), plot = gg, base_height = 2, base_width = 7)
