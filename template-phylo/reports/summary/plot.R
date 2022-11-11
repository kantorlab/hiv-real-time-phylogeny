library(assertthat)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(lubridate)
library(ggrepel)

args <- commandArgs(trailingOnly=TRUE)

study_file     <- args[1]
new_ids_file   <- args[2]
cluster_file   <- args[3]
treatment_file <- args[4]
dataset        <- args[5]
out_file       <- args[6]

# Data preprocessing

dataset <- substr(dataset, 1, 12)

new <- read.csv(new_ids_file) %>%
       mutate(New=TRUE) %>%
       select(StudyID, New)

treatment <- read_tsv(treatment_file) %>%
             mutate(Index=as.factor(ifelse(Dataset == dataset, "Current", "Previous"))) %>%
             select(StudyID, Index)

D <- left_join(read.csv(study_file), read.delim(cluster_file)) %>%
     left_join(new) %>%
     left_join(treatment) %>%
     mutate(HIVDxDate=mdy(str_extract(HIVDxDate, "\\d+/\\d+/\\d+")),
            HIVDxMonths=12*year(HIVDxDate)+month(HIVDxDate)+day(HIVDxDate)/31)

# Summarize recent clusters with three or more patients

summary <- D %>%
           filter(!is.na(ClusterID)) %>%
           group_by(ClusterID) %>%
           summarise(ClusterLast=max(HIVDxDate, na.rm=TRUE),
                     ClusterFirst=min(HIVDxDate, na.rm=TRUE),
                     ClusterNew=sum(New, na.rm=TRUE),
                     ClusterSize=n()) %>%
           arrange(ClusterLast, ClusterFirst, ClusterSize) %>%
           filter(ClusterLast > ymd("20170101")) %>%
           mutate(Order=row_number())
D_lumped <- D %>%
            filter(ClusterID %in% summary$ClusterID) %>%
            filter(!is.na(HIVDxDate)) %>%
            left_join(summary) %>%
            mutate(CaseType=factor(case_when(Index == "Current" ~ "Index Case - Current Month",
                                             Index == "Previous" ~ "Index Case - Previous Months",
                                             New ~ "New Sequence for a Previous Positive Case",
                                             TRUE ~ "Neither")))

# Re-center around 12 previous months

assert_that(all(!is.na(D_lumped$HIVDxMonths)))
max_months <- ceiling(max(D_lumped$HIVDxMonths))
print(max(D_lumped$HIVDxDate))
origin <- max_months - 11

# Produce the segment plot

rescale <- function(x) { ifelse(x > -1, x, -sqrt(-x)) }

D_lumped$HIVDxMonths <- rescale(D_lumped$HIVDxMonths - origin)

x_limits <- c(rescale(12*1980 + 1 - origin), 11)
print(x_limits)

x_labels <- c(1980, 2000, 2010, 2015, 2018, 2020, 2021)
x_breaks <- c(sapply(x_labels, function(months) { rescale(12*months + 1 - origin) } ))
x_labels <- c(x_labels, sapply(max_months - seq(12, 0), function(months) { paste0((months - 1) %/% 12, "-", str_pad((months - 1) %% 12 + 1, 2, "left", "0")) } ))
x_breaks <- c(x_breaks, seq(-1, 11))
print(x_labels)
print(x_breaks)

y_max <- max(D_lumped$Order)

gg_segment <- ggplot(D_lumped) +
              geom_line(aes(group=Order, x=HIVDxMonths, y=Order), col="grey60") +
              geom_point(aes(x=HIVDxMonths, y=Order, col=CaseType)) +
              geom_text_repel(
                   aes(x=HIVDxMonths, y=Order, label=StudyID),
                   data = filter(D_lumped, CaseType == "Index Case - Current Month"),
                   #direction = "y",
                   nudge_y = -10,
                   force = 10,
                   color = "grey60"
              ) +
              geom_bracket(xmin=-1, xmax=11, y.position=y_max+5, label="Previous 12 Months") +
              labs(title="Clusters since 2017") +
              scale_size_area(max_size=2) +
              theme_classic() +
              theme(axis.line=element_blank(),
                    axis.ticks.y=element_blank(),
                    axis.title=element_blank(),
                    axis.text.x=element_text(size=8, angle=90, vjust=0.5),
                    axis.text.y.right=element_text(size=7),
                    axis.text.y.left=element_blank(),
                    legend.justification=c(0, 1),
                    legend.position=c(0, 1),
                    panel.grid.major.x=element_line(colour="grey80", linetype="dashed", size=0.2)) +
              scale_x_continuous(limits=x_limits, breaks=x_breaks, labels=x_labels) +
              scale_colour_manual(name="This Month's Report", values=c("Index Case - Current Month"="red", "Index Case - Previous Months"="forestgreen", "New Sequence for a Previous Positive Case"="blue", "Neither"="grey60"))

pdf(out_file, width=13, height=13)
print(gg_segment)
dev.off()

