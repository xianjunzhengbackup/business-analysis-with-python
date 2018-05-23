# The R source code that drives the sensitivity analysis of the NPV from the
# risk model follows here.


# Read source data and function files. Modify the path names to match your
# directory structure and file names.
source("/Applications/R/RProjects/BizSimWithR/data/global_assumptions.R")
d.data <-
  read.csv("/Applications/R/RProjects/BizSimWithR/data/risk_assumptions.csv")
source("/Applications/R/RProjects/BizSimWithR/libraries/My_Functions.R")

# Slice the values from data frame d.data.
d.vals <- d.data[, 2:4]

sens.range <- c(0.1, 0.9)
len.d.vals <- length(d.vals[, 1])
len.sens.range <- length(sens.range)
npv.sens <- array(0, c(len.d.vals, len.sens.range))

# Assign values to variables using appropriate distributions.
p1.capex <- CalcBrownJohnson(0, d.vals[1, 1], d.vals[1, 2],
                             d.vals[1, 3], , kSampsize)
p1.dur <- round(CalcBrownJohnson(1, d.vals[2, 1], d.vals[2, 2],
                                 d.vals[2, 3], , kSampsize), 0)
p2.capex <- CalcBrownJohnson(0, d.vals[3, 1], d.vals[3, 2],
                             d.vals[3, 3], , kSampsize)
p2.dur <- round(CalcBrownJohnson(1, d.vals[4, 1], d.vals[4, 2],
                                 d.vals[4, 3], , kSampsize), 0)
maint.capex <- CalcBrownJohnson(0, d.vals[5, 1], d.vals[5, 2],
                                d.vals[5, 3], , kSampsize)
fixed.prod.cost <- CalcBrownJohnson(0, d.vals[6, 1], d.vals[6, 2],
                                    d.vals[6, 3], , kSampsize)
prod.cost.escal <- CalcBrownJohnson(, d.vals[7, 1], d.vals[7, 2],
                                    d.vals[7, 3], , kSampsize)
var.prod.cost <- CalcBrownJohnson(0, d.vals[8, 1], d.vals[8, 2],
                                  d.vals[8, 3], , kSampsize)
var.cost.redux <- CalcBrownJohnson(, d.vals[9, 1], d.vals[9, 2],
                                   d.vals[9, 3], , kSampsize)
gsa.rate <- CalcBrownJohnson(0, d.vals[10, 1], d.vals[10, 2],
                             d.vals[10, 3], 100, kSampsize)
time.to.peak.sales <- round(CalcBrownJohnson(1, d.vals[11, 1],
                                             d.vals[11, 2], d.vals[11, 3], , kSampsize), 0)
mkt.demand <- CalcBrownJohnson(0, d.vals[12, 1], d.vals[12, 2],
                               d.vals[12, 3], , kSampsize)
price <- CalcBrownJohnson(0, d.vals[13, 1], d.vals[13, 2],
                          d.vals[13, 3], , kSampsize)
rr.comes.to.market <- rbinom(kSampsize, 1, d.vals[14, 2])
rr.time.to.market <- round(CalcBrownJohnson(1, d.vals[15, 1],
                                            d.vals[15, 2], d.vals[15, 3], , kSampsize), 0)
early.market.share <- CalcBrownJohnson(0, d.vals[16, 1],
                                       d.vals[16, 2], d.vals[16, 3], 100, kSampsize)
late.market.share <- CalcBrownJohnson(0, d.vals[17, 1],
                                      d.vals[17, 2], d.vals[17, 3], 100, kSampsize)
price.redux <- CalcBrownJohnson(0, d.vals[18, 1], d.vals[18, 2],
                                d.vals[18, 3], , kSampsize)

d.vals.vect <- c(
  p1.capex,
  p1.dur,
  p2.capex,
  p2.dur,
  maint.capex,
  fixed.prod.cost,
  prod.cost.escal,
  var.prod.cost,
  var.cost.redux,
  gsa.rate,
  time.to.peak.sales,
  mkt.demand,
  price,
  rr.comes.to.market,
  rr.time.to.market,
  early.market.share,
  late.market.share,
  price.redux
)

d.vals.temp <- array(d.vals.vect, dim = c(kSampsize, len.d.vals))
d.vals.temp2 <- d.vals.temp
d.vals2 <- d.vals[,-2]

CalcBizSim = function(x) {
  # x is the data array that contains the pre-simulated samples for
  # each variable.

  p1.capex <- x[, 1]
  p1.dur <- x[, 2]
  p2.capex <- x[, 3]
  p2.dur <- x[, 4]
  maint.capex <- x[, 5]
  fixed.prod.cost <- x[, 6]
  prod.cost.escal <- x[, 7]
  var.prod.cost <- x[, 8]
  var.cost.redux <- x[, 9]
  gsa.rate <- x[, 10]
  time.to.peak.sales <- x[, 11]
  mkt.demand <- x[, 12]
  price <- x[, 13]
  rr.comes.to.market <- x[, 14]
  rr.time.to.market <- x[, 15]
  early.market.share <- x[, 16]
  late.market.share <- x[, 17]
  price.redux <- x[, 18]

  # CAPEX Module
  phase <- t(sapply(run, function(r)
    (year <= p1.dur[r]) * 1 +
      (year > p1.dur[r] & year <= (p1.dur[r] + p2.dur[r])) * 2 +
      (year > (p1.dur[r] + p2.dur[r])) * 3))

  capex <-
    t(sapply(run, function(r)
      (phase[r,] == 1) * p1.capex[r] / p1.dur[r] +
        (phase[r,] == 2) * p2.capex[r] / p2.dur[r] +
        (phase[r,] == 3) * maint.capex[r]))

  # Depreciation Module
  depr.matrix <-
    array(sapply(run, function(r)
      sapply(year, function(y)
        ifelse(
          y <= p1.dur[r] & year > 0,
          0,
          ifelse(
            y == (p1.dur[r] + 1) &
              year < y + kDeprPer & year >= y,
            p1.capex[r] / kDeprPer,
            ifelse((year >= y) &
                     (year < (y + kDeprPer)), capex[r, y - 1] / kDeprPer, 0)
          )
        ))),
      dim = c(kHorizon, kHorizon, kSampsize))

  depr <- t(sapply(run, function(r)
    sapply(year, function(y)
      sum(depr.matrix[y, , r]))))

  # Competition Module
  market.share <-
    (rr.comes.to.market == 1) * ((rr.time.to.market <= p1.dur) *
                                   early.market.share / 100 + (rr.time.to.market > p1.dur) *
                                   late.market.share / 100
    ) +
    (rr.comes.to.market == 0) * 1

  # Sales Module
  mkt.adoption <- t(sapply(run, function(r)
    market.share[r] *
      pmin(cumsum(phase[r,] > 1) / time.to.peak.sales[r], 1)))
  sales <-
    t(sapply(run, function(r)
      mkt.adoption[r,] * mkt.demand[r] *
        1000 * 2000))
  revenue <- t(sapply(run, function(r)
    sales[r,] * price[r] *
      (1 - rr.comes.to.market[r] * price.redux[r] / 100)))

  # OPEX Module
  fixed.cost <-
    t(sapply(run, function(r)
      (phase[r,] > 1) * fixed.prod.cost[r] *
        (1 + prod.cost.escal[r] / 100) ^ (year - p1.dur[r] - 1)))
  var.cost <- t(sapply(run, function(r)
    var.prod.cost[r] *
      (1 - var.cost.redux[r] / 100) ^ (year - p1.dur[r] - 1) * sales[r,]))
  gsa <- t(sapply(run, function(r)
    (gsa.rate[r] / 100) * revenue[r,]))
  opex <- fixed.cost + var.cost

  # Value
  gross.profit <- revenue - gsa
  op.profit.before.tax <- gross.profit - opex - depr
  tax <- op.profit.before.tax * kTaxRate / 100
  op.profit.after.tax <- op.profit.before.tax - tax
  cash.flow <- op.profit.after.tax + depr - capex
  cum.cash.flow <- t(sapply(run, function(r)
    cumsum(cash.flow[r,])))

  # Following the convention for when payments are counted as occurring
  # at the end of a time period.
  discount.factors <- 1 / (1 + kDiscountRate / 100) ^ year
  discounted.cash.flow <- t(sapply(run, function(r)
    cash.flow[r,] *
      discount.factors))
  npv <- sapply(run, function(r)
    sum(discounted.cash.flow[r,]))
  return(npv)
}

base.mean <- mean(CalcBizSim(d.vals.temp))

for (i in 1:len.d.vals) {
  for (k in 1:len.sens.range) {
    # For a given variable, replace its samples with a vector containing
    # each sensitivity end point.
    d.vals.temp2[, i] <-  rep(d.vals2[i, k], kSampsize)

    # Calculate the mean NPV by calling the CalcBizSim() function.
    mean.npv <- mean(CalcBizSim(d.vals.temp2))


    # Insert the resultant mean NPV into an array that catalogs the
    # variation in the mean NPV by each variables' sensitivity points.
    npv.sens[i, k] <- mean.npv
  }

  # Restore the current variable's last sensitivity point with its original
  # simulated samples.
  d.vals.temp2[, i] <- d.vals.temp[, i]
}

# Assign npv.sens to a data frame.
var.names <- d.data$variable
rownames(npv.sens) <- d.data$variable
colnames(npv.sens) <- sens.range

# Sets up the sensitivity array.
npv.sens.array <- array(0, c(len.d.vals, 2))
npv.sens.array[, 1] <- (npv.sens[, 1] - base.mean)
npv.sens.array[, 2] <- (npv.sens[, 2] - base.mean)
rownames(npv.sens.array) <- var.names
colnames(npv.sens.array) <- sens.range

# Calculates the rank order of the NPV sensitivity based on the
# absolute range caused by a given variable. The npv.sens.array
# is reindexed by this rank ordering for the bar plot.
npv.sens.rank <- order(abs(npv.sens.array[, 1] -
                             npv.sens.array[, 2]), decreasing = FALSE)
ranked.npv.sens.array <- npv.sens.array[npv.sens.rank,]
ranked.var.names <- var.names[npv.sens.rank]
rownames(ranked.npv.sens.array) <- ranked.var.names

# Plots the sensitivity array.
par(mai = c(1, 1.75, .5, .5))
barplot(
  t(ranked.npv.sens.array) / 1000,
  main = "NPV Sensitivity to
  Uncertainty Ranges",
  names.arg = ranked.var.names,
  col = "red",
  xlab = "NPV [$000]",
  beside = TRUE,
  horiz = TRUE,
  offset = base.mean / 1000,
  las = 1,
  space = c(-1, 1),
  cex.names = 1
)
