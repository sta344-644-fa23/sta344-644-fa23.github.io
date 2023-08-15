data(h02, package="fpp")
forecast::ggtsdisplay(h02, points=FALSE)

calc_rmse = function(m) {
  yardstick::rmse_vec(
    truth = unclass(m$x),
    estimate = unclass(m$fitted)
  )
}

# S Differencing

(m1 = forecast::Arima(
  h02, order = c(0,0,0), 
  seasonal = list(order=c(0,1,0), period=12)
))
# AICc=-452.11

forecast::ggtsdisplay(
  m1$residuals, points=FALSE
)

## S MA

(m2.1 = forecast::Arima(
  h02, order = c(0,0,0), 
  seasonal = list(order=c(0,1,1), period=12)
))

(m2.2 = forecast::Arima(
  h02, order = c(0,0,0), 
  seasonal = list(order=c(0,1,2), period=12)
))


calc_rmse(m1)
calc_rmse(m2.1)
calc_rmse(m2.2)

forecast::ggtsdisplay(
  m2.2$residuals, points=FALSE
)

## AR

(m3.1 = forecast::Arima(
  h02, order = c(1,0,0), 
  seasonal = list(order=c(0,1,2), period=12)
))

(m3.2 = forecast::Arima(
  h02, order = c(2,0,0), 
  seasonal = list(order=c(0,1,2), period=12)
))

(m3.3 = forecast::Arima(
  h02, order = c(3,0,0), 
  seasonal = list(order=c(0,1,2), period=12)
))

calc_rmse(m3.1)
calc_rmse(m3.2)
calc_rmse(m3.3)


forecast::ggtsdisplay(
  m3.3$residuals, points=FALSE
)

## Forecast

plot(forecast::forecast(m3.3))
lines(m3.3$fitted, col=adjustcolor('blue', alpha.f = 0.5))
