# Usage: Rscript 36_completion_components.R CDC.zip rebuilt_private.rds aggregate_output.csv
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3)
library(survey)
read_part <- function(name) {
  read.delim(unz(args[1], paste0("PartA03062025/", name, ".csv")), sep = "|", na.strings = "")
}
people <- read_part("PartAMain")
labs <- subset(read_part("AllCombinedLabs"), d_main_set == 1)
stopifnot(!anyDuplicated(people$studyid), !anyDuplicated(labs$studyid))
b <- subset(merge(people, labs, by = "studyid", all.x = TRUE), d_age_cat %in% 6:8)
b <- transform(b, foreign = as.integer(birth_country != "USA"), q = qft_pn %in% 0:1,
               t = tst_in_window %in% 1 & is.finite(tst_result) & tst_result >= 0)
n <- readRDS(args[2])$data
n <- transform(n, foreign = as.integer(birth == "NonUS"), q = LBXTBIN %in% 1:2,
               t = is.finite(TBDRUIND), history = TBQ070 %in% 2)
n <- subset(n, WTMEC2YR > 0)
nd <- svydesign(~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTMEC2YR, nest = TRUE, data = n)
nd <- subset(nd, age >= 25 & !(TBQ070 %in% 1) & !is.na(foreign))
bd <- svydesign(~d_site, weights = ~1, data = b)
bd <- subset(bd, !is.na(foreign))
compare <- function(design, endpoint, outcome, cohort) {
  design$variables$outcome <- as.numeric(outcome)
  model <- svyglm(outcome ~ foreign, design)
  est <- unname(coef(model))
  se <- sqrt(vcov(model)[2,2])
  width <- qt(.975, 17) * se
  data.frame(cohort, endpoint, us_pct = 100 * est[1], foreign_pct = 100 * sum(est),
             difference_pp = 100 * est[2], low = 100 * (est[2]-width), high = 100 * (est[2]+width))
}
components <- function(design, cohort) {
  d <- design$variables
  states <- list(paired = d$q & d$t, QFT_alone = d$q & !d$t, TST_alone = !d$q & d$t,
                 neither = !d$q & !d$t)
  tab <- do.call(rbind, lapply(names(states), function(s) compare(design, s, states[[s]], cohort)))
  stopifnot(abs(sum(tab$difference_pp)) < 1e-7)
  tab
}
out <- rbind(components(nd, "NHANES"), components(bd, "TBESC"),
             compare(nd, "negative_history_documented", nd$variables$history, "NHANES"),
             components(subset(nd, history), "NHANES documented history"))
write.csv(out, args[3], row.names = FALSE)
