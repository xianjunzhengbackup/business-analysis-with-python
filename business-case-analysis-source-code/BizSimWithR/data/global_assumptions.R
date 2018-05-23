# Global Model Assumptions
kHorizon <- 20 # yrs, length of time the model covers.
year <- 1:kHorizon # an index for temporal calculations.
kSampsize <- 1000 # the number of iterations in the Monte Carlo simulation.
run <- 1:kSampsize # the iteration index.
kTaxRate <- 38 # %
kDiscountRate <- 12 # %/yr, used for discounted cash flow calculations.
kDeprPer <- 7 # yrs, the depreciation schedule for the capital.
