library("ggplot2")
library("tidyverse")
library(RColorBrewer)
library(khroma)

setwd("/Users/Oscar/Desktop/Emerson_Lab_Work/SV-project")

ggplot(data.frame(x = c(0 , 200)), aes(x = x))+
#  scale_x_continuous(sec.axis = ~ . /9500) +
         stat_function(fun = dgamma, args = list(shape = 1.43, scale=38.98), aes(color="all_SV",)) +
  stat_function(fun = dgamma, args = list(shape = 0.159, scale=1378), aes(color="syn")) +
  stat_function(fun = dgamma, args = list(shape = 0.364, scale=1490), aes(color="nonsyn")) +
  labs(x = "Ne * s", y="frequency")