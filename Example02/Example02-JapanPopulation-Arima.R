# Step 0: Load libraries

library("forecast")
library("tseries")

# Step 1: Definition of time-series data from Japan

y_jap      <- c(124.3837, 125.9753, 127.0278, 127.8158, 128.1520, 127.9283, 127.1600)
ts_y_jap   <- ts(y_jap, start = 1992, frequency = 0.25)

# Step 2: Fit with auto.arima()

fit_jap <- auto.arima(ts_y_jap)

# Step 3: Plot forecast

forec_jap <- forecast(fit_jap, h = 15)
plot(forec_jap)