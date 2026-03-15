# Step 0: Load libraries

library("forecast")
library("tseries")

# Step 1: Definition of time-series data from Japan

y_jap      <- c(5675551259, 5842055733, 6007066690, 6171703001, 6337730340, 6503377776, 6671452017, 6844457659, 7021732156, 7201202481, 7381616240, 7558554525)
ts_y_jap   <- ts(y_jap, start = 1994, frequency = 0.5)

# Step 2: Fit with auto.arima()

fit_jap <- auto.arima(ts_y_jap)

# Step 3: Plot forecast

forec_jap <- forecast(fit_jap, h = 10)
plot(forec_jap, main="Forecast from ARIMA(1,2,0) for world's total population", xlab ="Time (Years)", ylab ="Population size")