# The R source code that drives the risk simulation follows here.

# Read source data and function files. Modify the path names to match your
# directory structure and file names.
source("/Applications/R/RProjects/BizSimWithR/data/global_assumptions.R")
d.data <-
  read.csv("/Applications/R/RProjects/BizSimWithR/data/risk_assumptions.csv")
source("/Applications/R/RProjects/BizSimWithR/libraries/My_Functions.R")

# Slice the values from data frame d.data.
d.vals <- d.data[, 2:4]

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
mean.npv <- mean(npv)

# Calculates the 80th percentile quantiles in the cash flow and
# cumulative cash flow.
q80 <- c(0.1, 0.5, 0.9)
cash.flow.q80 <-
  sapply(year, function(y)
    quantile(cash.flow[, y], q80))
cum.cash.flow.q80 <-
  sapply(year, function(y)
    quantile(cum.cash.flow[, y],
             q80))

# Plots the 80th percentile cash flow quantiles.
plot(
  0,
  type = "n",
  xlim = c(1, kHorizon),
  ylim = c(min(cash.flow.q80) / 1000,
           max(cash.flow.q80) / 1000),
  xlab = "Year",
  ylab = "[$000]",
  main = "Cash Flow",
  tck = 1
)
lines(
  year,
  cash.flow.q80[1, ] / 1000,
  type = "b",
  lty = 1,
  col = "blue",
  pch = 16
)
lines(
  year,
  cash.flow.q80[2,] / 1000,
  type = "b",
  lty = 1,
  col = "red",
  pch = 18
)
lines(
  year,
  cash.flow.q80[3,] / 1000,
  type = "b",
  lty = 1,
  col = "darkgreen",
  pch = 16
)
legend(
  "topleft",
  legend = q80,
  bg = "grey",
  pch = c(16, 18, 16),
  col = c("blue", "red", "dark green")
)

# Plots the 80th percentile cumulative cash flow quantiles.
plot(
  0,
  type = "n",
  xlim = c(1, kHorizon),
  ylim = c(min(cum.cash.flow.q80) / 1000,
           max(cum.cash.flow.q80) / 1000),
  xlab = "Year",
  ylab = "[$000]",
  main = "Cumulative Cash Flow",
  tck = 1
)
lines(
  year,
  cum.cash.flow.q80[1, ] / 1000,
  type = "b",
  lty = 1,
  col = "blue",
  pch = 16
)
lines(
  year,
  cum.cash.flow.q80[2, ] / 1000,
  type = "b",
  lty = 1,
  col = "red",
  pch = 18
)
lines(
  year,
  cum.cash.flow.q80[3, ] / 1000,
  type = "b",
  lty = 1,
  col = "darkgreen",
  pch = 16
)
legend(
  "topleft",
  legend = q80,
  bg = "grey",
  pch = c(16, 18, 16),
  col =
    c("blue", "red", "dark green")
)

# Calculates & plots the histogram of NPV.
breakpoints <-
  seq(min(npv), max(npv), abs(min(npv) - max(npv)) / 20)
hist(
  npv / 1000,
  freq = FALSE,
  breaks = breakpoints / 1000,
  main = "Histogram of NPV",
  xlab = "NPV [$000]",
  ylab = "Probability Density",
  col = "blue"
)

# Calculates the cumulative NPV probability chart and table.
cum.quantiles <- seq(0, 1, by = 0.05)
cum.npv.vals <-
  quantile(npv, cum.quantiles) #plot these values in chart
cum.npv.frame <- data.frame(cum.npv.vals) #table

# Plot the cumulative probability NPV curve.
plot(
  cum.npv.vals / 1000,
  cum.quantiles,
  main = "Cumulative Probability of NPV",
  xlab = "NPV [$000]",
  ylab = "Cumulative Probability",
  "b",
  tck = 1,
  col = "blue",
  pch = 16
)

# Pro Forma
# Create a data frame of the variables' mean values to be used in the
# pro forma.
pro.forma.vars <- array(
  c(
    sales,
    revenue,
    -gsa,
    gross.profit,-fixed.cost,
    -var.cost,
    -opex,
    -depr,
    op.profit.before.tax,
    -tax,
    op.profit.after.tax,
    depr,
    -capex,
    cash.flow
  ),
  dim = c(kSampsize, kHorizon, 14)
)

# Finds the annual mean of each pro forma element.
mean.pro.forma.vars <- array(0, c(14, kHorizon))

for (p in 1:14) {
  mean.pro.forma.vars[p,] <- sapply(year, function(y)
    mean(pro.forma.vars[, y, p]))
}

pro.forma <- data.frame(mean.pro.forma.vars)

# Assign text names to a vector. These will be the column headers of
# the data frame.
pro.forma.headers <-
  c(
    "Sales [lbs]",
    "Revenue",
    "GS&A",
    "Gross Profit",
    "Fixed Cost",
    "Variable Cost",
    "OPEX",
    "-Depreciation",
    "Operating Profit Before Tax",
    "Tax",
    "Operating Profit After Tax",
    "+Depreciation",
    "CAPEX",
    "Cash Flow"
  )

# Coerces the default column headers to be the headers we like.
colnames(pro.forma) <- year
rownames(pro.forma) <- pro.forma.headers

# Waterfall chart
# Extract the rows from the Pro Forma for the water fall chart.
waterfall.rows <- c(2, 3, 5, 6, 10, 13)
waterfall.headers <- pro.forma.headers[waterfall.rows]
wf.pro.forma <- pro.forma[waterfall.rows,]

# Find the present value of the extracted pro forma elements.
pv.wf.pro.forma <- rep(0, length(waterfall.rows))
for (w in 1:length(waterfall.rows)) {
  pv.wf.pro.forma[w] <- sum(wf.pro.forma[w,] * discount.factors)
}

cum.pv.wf.pro.forma1 <- cumsum(pv.wf.pro.forma)
cum.pv.wf.pro.forma2 <- c(0, cum.pv.wf.pro.forma1[1:(length(waterfall.rows) -
                                                       1)])
wf.high <- pmax(cum.pv.wf.pro.forma1, cum.pv.wf.pro.forma2)
wf.low <- pmin(cum.pv.wf.pro.forma1, cum.pv.wf.pro.forma2)

waterfall <- array(0, c(2, length(waterfall.rows)))
waterfall[1,] <- wf.high
waterfall[2,] <- wf.low
colnames(waterfall) <- waterfall.headers
rownames(waterfall) <- c("high", "low")

# Plot the waterfall.
boxplot(
  waterfall / 1000,
  data = waterfall / 1000,
  notch = FALSE,
  main = "Waterfall Chart",
  xlab = waterfall.headers,
  ylab = "$000",
  col = c("blue", rep("red", 5))
)
