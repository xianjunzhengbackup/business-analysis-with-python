# BrownJohnson distribution
CalcBrownJohnson = function(minlim=-Inf, p10, p50, p90, maxlim=Inf, n, samponly=TRUE) {
# This function simulates a distribution from three quantile
# estimates for the probability intervals of the predicted outcome.
# p10 = the 10th percentile estimate
# p50 = the 50th percentile estimate
# p90 = the 90th percentile estimate
# n = the number of runs used in the simulation
# Note, this function strictly requires that p10 < p50 < p90.
# The process of simulation is simple Monte Carlo.
# The returned result is, by default, a vector of values for X if
# samponly=TRUE, else a (n x 2) matrix is returned such that the first
# column contains the domain samples in X, and the second column
# contains the uniform variate samples from U.

  if (p10 == p50 && p50 == p90) {
		return(rep(p50, n))
	} else if (p10 >= p50 || p50 >= p90) {
    stop("Parameters not given in the correct order. Must be
    given as p10 < p50 < p90.")
  } else {
#Create a uniform variate sample space in the interval (0,1).
    U <- runif(n, 0, 1)

# Calculates the virtual tails of the distribution given the p10, p50, p90
# inputs. Truncates the tails at the upper and lower limiting constraints.
    p0 <- max(minlim, 2.5 * p10 - 1.5 * p50)
    p100 <- min(maxlim, 2.5 * p90 - 1.5 * p50)

#This next section finds the linear coefficients of the system of linear
# equations that describe the linear spline, using linear algebra...
# [C](A) = (X)
# (A) = [C]^-1 * (X)
# In this case, the elements of (C) are found using the values (0, 0.1,
# 0.5, 0.9, 1) at the end points of each spline segment. The elements
# of (X) correspond to the values of (p0, p10, p10, p50, p50, p90,
# p90, p100). Solving for this system of linear equations gives linear
# coefficients that transform values in U to intermediate values in X.
# Because there are four segments in the linear spline, and each
# segment contains two unknowns, a total of eight equations are
# required to solve the system.

# The spline knot values in the X domain.
    knot.vector <- c(p0, p10, p10, p50, p50, p90, p90, p100)

# The solutions to the eight equations at the knot points required to
# describe the linear system.
    coeff.vals <- c(0, 1, 0, 0, 0, 0, 0, 0, 0.1, 1, 0, 0, 0, 0, 0, 0,
    0, 0, 0.1, 1, 0, 0, 0, 0, 0, 0, 0.5, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0.5,
    1, 0, 0, 0, 0, 0, 0, 0.9, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0.9, 1, 0, 0, 0,
    0, 0, 0, 1, 1)

# The coefficient matrix created from the prior vector looks like the
# following matrix:
  #   [, 1]   [, 2]   [, 3]   [, 4]   [, 5]   [, 6]   [, 7]  [, 8]
# [1,] 0.0     1.0     0.0     0.0     0.0     0.0     0.0    0.0
# [2,] 1.0     1.0     0.0     0.0     0.0     0.0     0.0    0.0
# [3,] 0.0     0.0     1.0     1.0     0.0     0.0     0.0    0.0
# [4,] 0.0     0.0     0.5     1.0     0.0     0.0     0.0    0.0
# [5,] 0.0     0.0     0.0     0.0     0.5     1.0     0.0    0.0
# [6,] 0.0     0.0     0.0     0.0     0.9     1.0     0.0    0.0
# [7,] 0.0     0.0     0.0     0.0     0.0     0.0     0.9    1.0
# [8,] 0.0     0.0     0.0     0.0     0.0     0.0     1.0    1.0

    coeff.matrix <- t(matrix(coeff.vals, nrow=8, ncol=8))

# The inverse of the coefficient matrix.
    inv.coeff.matrix <- solve(coeff.matrix)

# The solution vector of the linear coefficients.
    sol.vect <- inv.coeff.matrix %*% knot.vector

#Builds the response by the piecewise linear sections
    X <- (U <= 0.1) * (sol.vect[1, 1] * U + sol.vect[2, 1]) +
      (U > 0.1 & U <= 0.5) * (sol.vect[3, 1] * U + sol.vect[4, 1]) +
      (U > 0.5 & U <= 0.9) * (sol.vect[5, 1] * U + sol.vect[6, 1]) +
      (U > 0.9 & U <= 1) * (sol.vect[7, 1] * U + sol.vect[8, 1])

    if (samponly == TRUE) {
      return(X)
    } else {
      X <- array( c(X, U), dim=c(n, 2))
      return(X)
    }

# Plots the array automatically after calculation.
# Comment out as necessary.
# plot(X[, 1], X[, 2], type="p", xlab="Outcomes", ylab="Cumulative
# Probability")
  }
}

# NPV Function
CalcNPV = function(series, time, dr, eotp = TRUE) {
  # series = the cash flow series
  # time = index over which series is allocated
  # dr = discount rate
  # eotp = end of time period calculation (default is TRUE)
  # or beginning of time period calculation (FALSE)

  this.NPV <- sum(series/(1 + dr)^(time - (1- eotp)))
}

# IRR Function
CalcIRR = function (series, time, irr0=0.1, tolerance=0.00001) {
  #Calculates the discount rate that produces a 0 NPV
  # of a series of end-of-year cash flows.
  # series = a vector of cash flows
  # time = the time array index
  # irr0 = the initial guess
  # tolerance = the error around 0 at which the goal seek stops

  irr <- c(irr0)
  npv <- CalcNPV(series, time, irr[1])
  i <- 1
  while ( abs(npv[i]) > tolerance ) {
    if (i ==1) {
      if ( npv[i] > 0 ) {
        (irr[i + 1] = irr[i] * 2)
      } else {
        (irr[i + 1] <- irr[i]/2)
      }
    } else {
      # uses Newton-Raphson method of convergence
      slope.npv <- (npv[i] - npv[i - 1])/
        (irr[i] - irr[i - 1])
      irr[i + 1] <- irr[i] - npv[i]/slope.npv
    }
    npv[i + 1] <- CalcNPV(series, time, irr[i + 1])
    i <- i + 1
  }
  irr <- irr[i]
}

# Equivalent Period Interest Rate Function
CalcEqPerIR = function(i, N) {
  # Calculates the equivalent interest rate for subperiod of a stated period
  # interest rate
  # i = the stated period interest rate
  # N = then number of subperiods within the stated period interest rate

  eq.per.ir = ((1+i)^(1/N))-1
}

# Period Payment Function
CalcPerPayment = function(i, L, N) {
  # Calculates the period payment due on a loan amount, L, to be paid back at
  # the number of payment periods specified, N, at the period interest
  # rate i.

  per.payment = L*i/(1-(1+i)^(-N))
}

# Period Principal Payment Function
CalcPerPrinPayment = function(i, L, M, N) {
  # Returns the period N principal payments due on a loan L obtained at N=0.
  # The loan has a life of M periods paid back at a period interest rate i.
  # N should be a sequence from 1 to M.

  per.pr.in.payment = L*i*(1/(1-(1+i)^(-M))-1)*(1+i)^(N-1)
}

# Period Interest Payment Function
CalcPerInPayment = function(i, L, M, N) {
  # Returns the period N interest payments due on a loan L obtained at N=0.
  # The loan has a life of M periods paid back at a period interest rate i.
  # N should be a sequence from 1 to M.

  per.in.payment = L*i*((1+i)^(N-1)-((1+i)^(N-1)-1)/(1-(1+i)^(-M)))
}

# Period Remaining Balance Function
CalcPerRemBal = function(i, L, M, N) {
  # Returns the period N remaining balance due on a loan L obtained at N=0.
  # The loan has a life of M periods paid back at a period interest rate i.
  # N should be a sequence from 1 to M.

  per.rem.bal = L*((1+i)^N-((1+i)^N-1)/(1-(1+i)^(-M)))
}

# Future Worth of an Annuity Function
CalcFutureWorthAnnuity = function(A,i, N) {
  # The future worth of an annuity for a given period of time is the sum at
  # the end of that period of the future worths of all payments at a given
  # rate of interest each compounding period.
  # A = the amount of the annuity
  # i = the period interest rate
  # N = the period at the end of the compounding

  future.worth.annuity = (A*(1+i)^N-1)/i
}

# Present Worth of an Annuity Function
CalcPrsWorthAnnuity = function(A,i, N) {
  # The present worth of an annuity for a given period of time is the sum at
  # the start of that period of the present worths of each payment in the
  # series, at a given rate of interest each compounding period.
  # interest each compounding period.
  # A = the amount of the annuity
  # i = the period interest rate
  # N = the period at the end of the compounding

  prs.worth.annuity = A*((1+i)^N-1)/(i*(1+i)^N)
}

# Annuity from a Present Amount Function
CalcAnnuityPrsAmount = function(P,i, N) {
  # The annuity from a present amount is the period amount which can be
  # withdrawn for a definite period of time from a present sum of money
  # at a given rate of interest each compounding period.
  # P = the present value
  # i = the period interest rate
  # N = the period at which the present amount is depleted

  annuity.present.amount = P*i*(1+i)^N/((1+i)^N-1)
}

# Weighted Average Cost of Capital Function
CalcWACC = function(Rt, Ri, Rf, Rcs, Beta, Debt, Equity) {
  # Calculates the annual nominal (adjusted for inflation) weighted average
  # cost of capital of an investment.
  # Rt = Corporate tax rate
  # Ri = Interest rate on debt
  # Rf = Risk free bond rate. Typically 6.4% on government long-term bonds.
  # Rcs = Historical return on common stock.  Typically 7% (1926-1995).
  # Beta = Equity beta
  # Debt = Outstanding leverage
  # Equity = Market value of equity

  wacc = ((1-Rt)*Ri*Debt+(Rf+(Beta*Rcs))*Equity)/(Debt+Equity)
}

# Black-Scholes Option Value Function
CalcOptionVal = function(P,S,irf,vol,t,sel) {
  # P = current price of the asset
  # S = strike price of the option
  # irf = risk free interest rate
  # vol = annual volatility
  # t = time to expiration
  # sel = a flag to indicate whether the option considered is a put
  # or a call.
  # call = 1; put = 0

  sel = 1-selector
  d1 = (log(P/S) + (irf + (vol^2)/2) * t) / (vol * sqrt(t))
  d2 = d1 - (vol * sqrt(t))
  N1 = pnorm(d1 * (-1)^(sel), 0 , 1)
  N2 = pnorm(d2 * (-1)^(sel), 0 , 1)

  option.val = (-1)^(sel) * (P * N1 - S * exp(-irf * t) * N2)
}
