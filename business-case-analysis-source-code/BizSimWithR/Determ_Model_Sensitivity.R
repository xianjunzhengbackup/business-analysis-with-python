# The R file that drives the deterministic sensitivity analysis of the NPV follows here.

# Read source data and function files. Modify the path names to match your
# directory structure and file names.
source("/Applications/R/RProjects/BizSimWithR/data/global_assumptions.R")
d.data <-
  read.csv("/Applications/R/RProjects/BizSimWithR/data/risk_assumptions.csv")
source("/Applications/R/RProjects/BizSimWithR/libraries/My_Functions.R")

# Slice the values from data frame d.data.
d.vals <- d.data$p50[1:13]

d.vals.sens.points <- d.data[1:13, 2:4]
sens.point <- 1:3
len.d.vals <- length(d.vals)
len.sens.range <- length(sens.point)
npv.sens <- array(0, c(len.d.vals, len.sens.range))

for (i in 1:len.d.vals) {
  for (k in 1:len.sens.range) {
    d.vals[i] <-  d.vals.sens.points[i, k]

    # Assign values to variables.
    p1.capex <- d.vals[1]
    p1.dur <- d.vals[2]
    p2.capex <- d.vals[3]
    p2.dur <- d.vals[4]
    maint.capex <- d.vals[5]
    fixed.prod.cost <- d.vals[6]
    prod.cost.escal <- d.vals[7]
    var.prod.cost <- d.vals[8]
    var.cost.redux <- d.vals[9]
    gsa.rate <- d.vals[10]
    time.to.peak.sales <- d.vals[11]
    mkt.demand <- d.vals[12]
    price <- d.vals[13]

    # CAPEX Module
    phase <- (year <= p1.dur) * 1 +
      (year > p1.dur & year <= (p1.dur + p2.dur)) * 2 +
      (year > (p1.dur + p2.dur)) * 3

    capex <- (phase == 1) * p1.capex / p1.dur +
      (phase == 2) * p2.capex / p2.dur +
      (phase == 3) * maint.capex

    # Depreciation Module
    depr.matrix <-
      t(sapply(year, function(y)
        ifelse(
          y <= p1.dur & year > 0,
          0,
          ifelse(
            y == (p1.dur + 1) & year < y + kDeprPer & year >= y,
            p1.capex / kDeprPer,
            ifelse((year >= y) & (year < (y + kDeprPer)),
                   capex[y - 1] / kDeprPer, 0)
          )
        )))
    depr <- colSums(depr.matrix)

    # Sales Module
    mkt.adoption <- pmin(cumsum(phase > 1) / time.to.peak.sales, 1)
    sales <- mkt.adoption * mkt.demand * 1000 * 2000
    revenue <- sales * price

    # OPEX Module
    fixed.cost <- (phase > 1) * fixed.prod.cost *
      (1 + prod.cost.escal / 100) ^ (year - p1.dur - 1)

    var.cost <- var.prod.cost * (1 - var.cost.redux / 100) ^
      (year - p1.dur - 1) * sales

    gsa <- (gsa.rate / 100) * revenue
    opex <- fixed.cost + var.cost

    # Value
    gross.profit <- revenue - gsa
    op.profit.before.tax <- gross.profit - opex - depr
    tax <- op.profit.before.tax * kTaxRate / 100
    op.profit.after.tax <- op.profit.before.tax - tax
    cash.flow <- op.profit.after.tax + depr - capex
    cum.cash.flow <- cumsum(cash.flow)

    discount.factors <- 1 / (1 + kDiscountRate / 100) ^ year
    # Following the convention for when payments are counted as occurring
    # at the end of a time period.
    discounted.cash.flow <- cash.flow * discount.factors
    npv <- sum(discounted.cash.flow)
    npv.sens[i, k] <- npv
  }
  d.vals[i] <- d.data$p50[i]
}

# Assign npv.sens to a data frame.
var.names <- d.data[1:13, 1]
sens.point.names <- c("p10", "p50", "p90")
rownames(npv.sens) <- var.names
colnames(npv.sens) <- sens.point.names

# Sets up the sensitivity array.
npv.sens.array <- array(0, c(len.d.vals, 2))
npv.sens.array[, 1] <- (npv.sens[, 1] - npv.sens[, 2])
npv.sens.array[, 2] <- (npv.sens[, 3] - npv.sens[, 2])
rownames(npv.sens.array) <- var.names
colnames(npv.sens.array) <- sens.point.names[-2]

# Calculates the rank order of the NPV sensitivity based on the
# absolute range caused by a given variable. The npv.sens.array
# is reindexed by this rank ordering for the bar plot.
npv.sens.rank <-
  order(abs(npv.sens.array[, 1] - npv.sens.array[, 2]),
        decreasing = FALSE)

ranked.npv.sens.array <- npv.sens.array[npv.sens.rank,]
ranked.var.names <- var.names[npv.sens.rank]
rownames(ranked.npv.sens.array) <- ranked.var.names

# Plots the sensitivity array.
par(mai = c(1, 1.75, .5, .5))
barplot(
  t(ranked.npv.sens.array) / 1000,
  main = "Deterministic NPV
  Sensitivity",
  names.arg = ranked.var.names,
  col = "light blue",
  xlab = "NPV [$000]",
  beside = TRUE,
  horiz = TRUE,
  offset = npv.sens[, 2] / 1000,
  las = 1,
  space = c(-1, 1),
  cex.names = 1,
  tck = 1
)
