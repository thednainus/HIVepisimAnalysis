#source observed data
source(system.file("data/incidence_HIVdiagnosis.R", package = "HIVepisimAnalysis"))
incidence <- readRDS(system.file("data/ECDC_incidence_model_22Oct2021.RDS",
                                 package = "HIVepisimAnalysis"))

library(ggplot2)
library(reshape2)

#plot only the data in which migration rate was set to 500 migrants per year


param_1067_500 <- readRDS("Results_paper/all_diag_m_and_q_1067_500migrants.RDS")
param_1067_500["param"] <- "1067"
param_1067_500["migrant"] <- "500"

param_2348_500 <- readRDS("Results_paper/all_diag_m_and_q_2348_500migrants.RDS")
param_2348_500["param"] <- "2348"
param_2348_500["migrant"] <- "500"


all_diag <- rbind(param_1067_500[6:41,],
                  param_2348_500[6:41,])
all_diag["param_migrant"] <- paste(all_diag$param, all_diag$migrant, sep = "_")

all_diag$year <- as.character(all_diag$year)
all_diag$year <- as.numeric(all_diag$year)
all_diag$param_migrant <- as.factor(all_diag$param_migrant)

all_diag <- all_diag[,c(1:4,7)]

names(all_diag)[2:4] <- c("lower", "median", "upper")
all_diagm <- melt(all_diag, id.vars = c("year", "lower", "upper", "param_migrant"))


quartz()
names(incidenceDiag)[1] <- "year"
incidenceDiag$year <- as.numeric(incidenceDiag$year)
incidenceDiag <- incidenceDiag[6:41,]
ggplot(all_diagm, aes(x=year)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = param_migrant), alpha=0.20) +
  geom_line(aes(y = value, linetype = variable, colour = param_migrant), linetype = 2) +
  theme_bw() + ylab("Incidence of diagnosis") +
  scale_fill_manual(values=c("#c9222a", "#222ac9"), guide = "none") +
  scale_color_manual(values=c("#c9222a", "#222ac9", "black"),
                     breaks=c("1067_500", "2348_500", "San Diego data"),
                     labels=c("Parameters 1", "Parameters 2", "San Diego data")) +
  theme(text = element_text(size=20), legend.position = "bottom") +
  theme(legend.title=element_blank())

quartz()
names(incidenceDiag)[1] <- "year"
incidenceDiag$year <- as.numeric(incidenceDiag$year)

ggplot(all_diagm, aes(x=year)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = param_migrant), alpha=0.20) +
  geom_line(aes(y = value, linetype = variable, colour = param_migrant), linetype = 2) +
  geom_line(data = incidenceDiag[6:41,], aes(y = frequency, colour = "San Diego data"),
            size = 0.8) +
  theme_bw() + ylab("Incidence of diagnosis") +
  scale_fill_manual(values=c("#c9222a", "#222ac9"), guide = "none") +
  scale_color_manual(values=c("#c9222a", "#222ac9", "black"),
                     breaks=c("1067_500", "2348_500", "San Diego data"),
                     labels=c("Parameters 1", "Parameters 2", "San Diego data")) +
  theme(text = element_text(size=20), legend.position = "bottom") +
  theme(legend.title=element_blank())

ggsave("incidence_diagnosis_50migrants.pdf", useDingbats=FALSE, width=12, height=9)


#read data from other cities
NYC <- read.csv("../../Data_other_USA_areas/NYC.csv")
NYC <- NYC[c(1:36),]
names(NYC)[1] <- "year"
NYC$year <- as.numeric(NYC$year)
NYC$HIV.Dx.per.year <- as.numeric(NYC$HIV.Dx.per.year)

Seattle <- read.csv("../../Data_other_USA_areas/Seattle.csv")
Seattle <- Seattle[c(1:36),c(1:2)]
names(Seattle)[1] <- "year"
Seattle$year <- as.numeric(Seattle$year)
Seattle$HIV.Dx.per.year <- as.numeric(Seattle$HIV.Dx.per.year)

Miami <- read.csv("../../Data_other_USA_areas/Miami.csv")
Miami <- Miami[c(1:36),c(1:2)]
names(Miami)[1] <- "year"
Miami$year <- as.numeric(Miami$year)
Miami$HIV.Dx.per.year <- as.numeric(Miami$HIV.Dx.per.year)

Chicago <- read.csv("../../Data_other_USA_areas/Chicago.csv")
Chicago <- Chicago[c(1:36),c(1:2)]
names(Chicago)[1] <- "year"
Chicago$year <- as.numeric(Chicago$year)
Chicago$HIV.Dx.per.year <- as.numeric(Chicago$HIV.Dx.per.year)


quartz()
names(incidenceDiag)[1] <- "year"
incidenceDiag$year <- as.numeric(incidenceDiag$year)

ggplot(all_diagm, aes(x=year)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = param_migrant), alpha=0.20) +
  geom_line(aes(y = value, linetype = variable, colour = param_migrant), linetype = 2) +
  geom_line(data = incidenceDiag[6:41,], aes(y = frequency, colour = "San Diego"),
            size = 0.8) +
  geom_line(data = NYC, aes(y = HIV.Dx.per.year, colour = "NYC"),
            size = 0.8) +
  geom_line(data = Seattle, aes(y = HIV.Dx.per.year, colour = "Seattle"),
            size = 0.8) +
  geom_line(data = Chicago, aes(y = HIV.Dx.per.year, colour = "Chicago"),
            size = 0.8) +
  geom_line(data = Miami, aes(y = HIV.Dx.per.year, colour = "Miami"),
            size = 0.8) +
  theme_bw() + ylab("Incidence of diagnosis") +
  scale_fill_manual(values=c("#c9222a", "#222ac9"), guide = "none") +
  scale_color_manual(values=c("#c9222a", "#222ac9",
                              "black", "blue", "magenta", "gray", "orange"),
                     breaks=c("1067_500", "2348_500",
                              "San Diego", "NYC", "Seattle", "Chicago", "Miami"),
                     labels=c("Parameters 1", "Parameters 2",
                              "San Diego", "NYC", "Seattle", "Chicago", "Miami")) +
  theme(text = element_text(size=20), legend.position = "bottom") +
  theme(legend.title=element_blank())


