# OBTENCION Y ANALISIS EXPLORATORIO DE LOS DATOS
getwd()
setwd("D:/Yo/Universidad/TECI/SERIES TEMPORALES")
ruta_archivo <- "raintemp.csv"
datos <- read.csv(ruta_archivo)
head(datos)

datos_meses <- datos[datos$Type.of.period=='Monthly',]
head(datos_meses)

year_mon <- paste(datos_meses$Year, datos_meses$Period, sep = "-")
year_mon

### RESERVAMOS UN CONJUNTO DE TEST Y OTRO DE TRAIN
Temp <- datos_meses$Avg.temp.in.centigrade.
Temp_train <- Temp[1:(length(Temp)-12)]
Temp_test <- Temp[(length(Temp)-11):length(Temp)]
head(Temp)

## REPRESENTACI�N SERIE
plot(seq_along(Temp_train), Temp_train, type="l", col="blue", xlab="�ndice", ylab="Temperatura Media", main="Serie Temporal")

## IMPORTAMOS LIBRERIAS NECESARIAS
library(forecast)
library(lmtest)
library(timsac)
library(haven)
library(descomponer)
library(tsoutliers)
library(expsmooth)
library(fma)
library(tseries)

### DESCOMPOSICI�N SERIE ENTERA
datos <- ts(Temp_train,frequency = 12)
descomp <- decompose(datos)
plot(descomp)

#PASAMOS A KELVIN PARA HACER BOX-COX
Kelvin <- Temp_train +273.15
datos_kelvin<- ts(Kelvin,frequency = 12)
(lambda = forecast::BoxCox.lambda(datos_kelvin))
datos_BoxCox = forecast::BoxCox(datos_kelvin,lambda)
plot(datos_BoxCox)
descomp_BoxCox <- decompose(datos_BoxCox)
plot(descomp_BoxCox)
# NO HACE NADA LA TRANSFORMACION

#Detectar Outliers
outliers.datos = tsoutliers::tso(datos,types="AO",maxit.iloop=10)
outliers.datos
plot(outliers.datos)
#No hay outliers


### DECIDIMOS RECORTAR EL CONJUNTO A PARTIR DEL A�O 4
Temp_train_acortado <- Temp[49:(length(Temp)-12)]
head(Temp_train_acortado)

#Grafico
plot(seq_along(Temp_train_acortado), Temp_train_acortado, type="l", col="blue", xlab="�ndice", ylab="Temperatura Media", main="Serie Temporal")

#Descomposicion
datos_train_acortados <- ts(Temp_train_acortado,frequency = 12)
descomp <- decompose(datos_train_acortados)
plot(descomp)


#OBSERVAMOS COMPONENTE ESTACIONAL CON PERIODO 12
#Aplicamos diferenciaci�n estacional
(ndifes=nsdiffs(datos_train_acortados))
dif_es=diff(datos_train_acortados, lag=12, differences = ndifes)
plot(dif_es)
descomp_diffs <- decompose(dif_es)
plot(descomp_diffs)

#No hay componente regular
(ndif=ndiffs(dif_es))
dif_es_sim=diff(dif_es, lag=1, differences = ndif)
plot(dif_es_sim)
descomp_diffs_sim <- decompose(dif_es_sim)
plot(descomp_diffs_sim)


#Calculamos y mostramos el periodograma de la serie
#transformada y desestacionalizada
periodograma(dif_es)
gperiodograma(dif_es)

#Calculamos y mostramos el periodograma de la serie original
periodograma(datos_train_acortados)
gperiodograma(datos_train_acortados)

#Realizamos el test de Dickey-Fuller para la estacionariedad
resultado_prueba <- adf.test(dif_es)
print(resultado_prueba)


#VEAMOS LOS PAR�METROS DEL MODELO.

#Calculamos la FAS
acf(dif_es,lag.max = 50)

#Calculamos la FAP
pacf(dif_es,lag.max = 50)

#A partir de ahora suponemos que no hay interveciones 
#Tambien existe auto.arima que te prueba distintas variables.
auto.arima(datos_train_acortados, allowdrift=F,trace=T)

#Ajustamos los modelos m�s adecuados:
fitARIMA <- arima(datos_train_acortados, order=c(0,0,1), seasonal = list(order = c(1,1,0),
                                                                         period = 12),method="ML")
coeftest(fitARIMA) #Significaci�n par�metros
confint(fitARIMA) #Intervalos de confianza cada par�metro
summary(fitARIMA) #Coeficientes del modelo y m�tricas de error


fitARIMA2 <- arima(datos_train_acortados, order=c(0,0,1), seasonal = list(order = c(0,1,1),
                                                                         period = 12),method="ML")
coeftest(fitARIMA2) 
confint(fitARIMA2) 
summary(fitARIMA2)                                                                      


fitARIMA3 <- arima(datos_train_acortados, order=c(1,0,0), seasonal = list(order = c(1,1,0),
                                                                          period = 12),method="ML")
coeftest(fitARIMA3) 
confint(fitARIMA3) 
summary(fitARIMA3)                                                                        


fitARIMA4 <- arima(datos_train_acortados, order=c(1,0,0), seasonal = list(order = c(0,1,1),
                                                                          period = 12),method="ML")
coeftest(fitARIMA4) 
confint(fitARIMA4) 
summary(fitARIMA4)                                                                        


#PRUEBAS DE AMPLIACION DE MODELO
fitARIMA_amp <- arima(datos_train_acortados, order=c(0,0,2), seasonal = list(order = c(0,1,1),
                                                                          period = 12),method="ML")
coeftest(fitARIMA_amp) 
confint(fitARIMA_amp) 
summary(fitARIMA_amp) 


# Prueba de Independencia: Ljung-Box
checkresiduals(fitARIMA2) # S� cumple independencia (p-valor >0.05)

#Homocedasticidad:
plot(fitARIMA2$residuals)

#Normalidad->
#Grafica:
qqnorm(fitARIMA2$residuals)
qqline(fitARIMA2$residuals)

#Tests:
#shapiro.test(fitARIMA2$residuals)
#install.packages("nortest")
library(nortest)
lillie.test(fitARIMA2$residuals) #p-valor >0.05 

#Predicciones->

predict(fitARIMA2,n.ahead=12) #se usa m�s predict

futurVal <- forecast(fitARIMA2, h=36,level = c(95)) #alpha = 0.05 (ajuste de predicci�n al 95%)
plot(futurVal)


# Concatenar los conjuntos de datos
Temp_total <- c(Temp_train_acortado, Temp_test)

Pred_total <- c(Temp_train_acortado,futurVal$mean)

# Crear el gr�fico con los datos originales

plot(seq_along(Temp_test), Temp_test, type="l", col="blue", xlab="�ndice", ylab="Temperatura Media", main="Serie Temporal")
lines(seq_along(futurVal$mean), futurVal$mean, col="red", lty=2) 
legend("topright", legend=c("Real", "Predicci�n"), col=c("blue", "red"), lty=c(1, 2), cex=0.7)

# Calcular el error cuadr�tico medio (RMSE)
rmse <- sqrt(mean((Temp_test - futurVal$mean)^2))

# Imprimir el RMSE
cat("Error Cuadr�tico Medio (RMSE):", rmse, "\n")
