# Autocallable Certificate Pricing

**Financial Engineering** · Politecnico di Milano · Anna Belli · May 2026

---

## Overview

This project prices an **autocallable certificate** issued by Bank XX on the **EURO STOXX 50** index (100M EUR notional). The product features a conditional digital coupon and an early redemption clause: if the index closes below a strike of K = 3200 at the Year 1 reset date, a 6% coupon is paid and the note redeems early; otherwise, the contract runs to maturity and a 2% coupon is paid.

To hedge the structured payoff, Bank XX enters a swap with an Investment Bank. The goal is to determine the fair **upfront payment X%** that makes the NPV of the swap equal to zero at inception.

The valuation is set as of **15 February 2008**, using EUR market data (deposits, futures, and swaps) to bootstrap the discount curve.

The valuation is carried out under three frameworks:
- **NIG model** (Normal Inverse Gaussian) — via Lewis analytical formula and Monte Carlo simulation
- **VG model** (Variance Gamma) — via Monte Carlo simulation
- **Black-76 model** — plain and smile-adjusted, as a benchmark

---

## Repository Structure

```
Certificate-Pricing/
├── README.md
├── report/
│   └── Certificate_Pricing_report.pdf
├── data/
│   └── eurostoxx_Poli.mat           # EURO STOXX 50 implied vol surface
├── utilities_bootstrap/             # Market data and curve bootstrap helpers
│   ├── MktData_CurveBootstrap.xls   # EUR market data (deposits, futures, swaps)
│   ├── bootstrap.m                  # Discount curve bootstrap (depos + futures + swaps)
│   ├── readExcelData_mac.m          # Reads bid/ask rates and dates from Excel for Mac
│   ├── readExcelData_windows.m      # Reads bid/ask rates and dates from Excel for Windows
│   ├── fromdatetodiscount.m         # Interpolates zero rates → discount factors
│   └── ConvertDates.m               # Adjusts dates to business days (Modified Following)
├── utilities_exercise/              # Pricing functions for project
│   ├── calibrate.m                  # NIG/VG calibration via FFT + fmincon
│   ├── conditionOnEta.m             # Nonlinear constraint on eta parameter
│   ├── FFT_obj.m                    # European call pricing via Lewis + FFT
│   ├── plotting_distribution.m      # Subordinator density comparison (IG vs Gamma)
│   ├── UpfrontPricingMC.m           # Monte Carlo upfront pricing (NIG & VG)
│   ├── UpfrontPricingLEWIS.m        # Analytical upfront pricing via Lewis formula
│   ├── UpfrontPricingGamma.m        # Monte Carlo upfront pricing (VG-specific)
│   └── UpfrontPricingBS.m           # Black-76 upfront pricing with smile correction
└── RunExercise.m            # Main script
```

---

## Methods

### 1. Curve Bootstrap
The valuation date is **15 February 2008**. The discount curve is bootstrapped from EUR market data using deposits (short end), futures (medium term), and swaps (long end) via `bootstrap.m`. The input data are stored in `MktData_CurveBootstrap.xls`. Discount factors at any target date are then obtained via linear interpolation of zero rates (`fromdatetodiscount.m`).

### 2. Calibration
NIG and VG parameters (σ, κ, η) are calibrated to the EURO STOXX 50 implied volatility surface via FFT-based option pricing (`FFT_obj.m`) and `fmincon` least-squares minimisation (`calibrate.m`).

| Model | σ | κ | η |
|-------|--------|--------|---------|
| NIG | 0.1742 | 0.4526 | 10.2141 |
| VG | 0.2115 | 0.6674 | 4.2850 |

### 3. Valuation
The key quantity is the risk-neutral digital put probability P(S_{t1} ≤ K), which drives both the coupon leg and the floating leg NPV.

**Lewis formula (NIG):** the digital put probability is computed analytically via adaptive quadrature (`quadgk`), providing a closed-form benchmark.

**Monte Carlo (NIG & VG):** 10⁷ paths are simulated using the Inverse Gaussian (NIG) or Gamma (VG) subordinator. Boolean path indicators handle the early redemption logic. The 3-year extension simulates two consecutive increments to capture the joint distribution of (S_{t1}, S_{t2}).

**Black-76:** the digital probability is N(−d₂), with an optional smile correction via central finite difference on the vol surface (slope impact term).

---

## Results

### 2-year maturity

| Method | Model | Digital Put Prob. | Upfront X% |
|--------|-------|------------------|-----------|
| Lewis formula | NIG | 24.71% | 5.8595% |
| Monte Carlo | NIG | 24.73% | 5.8583% |
| Monte Carlo | VG | 25.06% | 5.8296% |
| Black-76 plain | — | — | 5.2142% |
| Black-76 smile-adjusted | — | — | 5.8275% |

The smile correction on Black-76 accounts for ~61 bps and brings the result in line with the Lévy models, confirming the importance of the volatility skew for digital payoffs.

### 3-year maturity (Monte Carlo only)

| Model | Upfront X% |
|-------|-----------|
| NIG | 8.1761% |
| VG | 8.1291% |

The 3-year extension rules out Black-76, which cannot capture the path-dependency across two reset dates.

---

## Requirements

- MATLAB R2021b or later
- Financial Toolbox (`blkprice`, `yearfrac`)
- Optimization Toolbox (`fmincon`)
- Statistics and Machine Learning Toolbox (`random('InverseGaussian', ...)`)

---

## How to Run

1. Open MATLAB.
2. Run the main script:
```matlab
run('RunExercise.m')

The script loads market data, bootstraps the discount curve, calibrates the models, and prints all upfront results to the console.

---

## References

- Baviera R., *Financial Engineering* lecture notes, Politecnico di Milano  
- Lewis A. (2001), *A Simple Option Formula for General Jump-Diffusion and Other Exponential Lévy Processes*
