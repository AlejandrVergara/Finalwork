library(tidyverse)
library(caret)
library(class)
library(gmodels)
library(psych)
library(DynamicCancerDriverKM)

view(DynamicCancerDriverKM::BRCA_normal)
datanormal <-(DynamicCancerDriverKM::BRCA_normal)
dataPt <-(DynamicCancerDriverKM::BRCA_PT)


# Calcula el porcentaje de valores menores a 10 en cada columna
porcentaje_menor_10 <- final_data %>%
  summarise_all(~ mean(. <100, na.rm = TRUE))

# Filtra las columnas donde el 80% o más de los valores son menores a 10
columnas_a_eliminar <- names(porcentaje_menor_10[, porcentaje_menor_10 >= 0.8])

# Elimina las columnas seleccionadas
final_data_filtrado <- final_data %>%
  select(-one_of(columnas_a_eliminar))




