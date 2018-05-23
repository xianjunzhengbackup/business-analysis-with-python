# The R source code that drives the initial deterministic simulation follows here.

# Read source data and function files. Modify the path names to match your
# directory structure and file names.
source("/Applications/R/RProjects/BizSimWithR/data/global_assumptions.R")
d.data <-
  read.csv("/Applications/R/RProjects/BizSimWithR/data/risk_assumptions.csv")
source("/Applications/R/RProjects/BizSimWithR/libraries/My_Functions.R")

# Slice the values from data frame d.data.
d.vals <- d.data$p50[1:13]

# Assign p50 values to variables.
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
    ))
  )
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

# Pro Forma
# Create an array of the variables to be used in the pro forma.
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
  dim = c(kHorizon, 14)
)

# Create a data frame of the pro forma array.
pro.forma <- data.frame(pro.forma.vars)

# Assign text names to a vector. These will be the column headers of the
# data frame.
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
colnames(pro.forma) <- pro.forma.headers
rownames(pro.forma) <- year

# Transposes the data frame so that columns become rows.
pro.forma = t(pro.forma)
