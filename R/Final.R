library(tidyverse)
library(caret)
library(class)
library(gmodels)
library(psych)
library(DynamicCancerDriverKM)

datanormal <-(DynamicCancerDriverKM::BRCA_normal)
dataPt <-(DynamicCancerDriverKM::BRCA_PT)
final_data <- bind_rows(datanormal, dataPt)


# Calcula el porcentaje de valores menores a 10 en cada columna
porcentaje_menor_10 <- final_data %>%
  summarise_all(~ mean(. <100, na.rm = TRUE))

# Filtra las columnas donde el 80% o más de los valores son menores a 10
columnas_a_eliminar <- names(porcentaje_menor_10[, porcentaje_menor_10 >= 0.8])

# Elimina las columnas seleccionadas
final_data_filtrado <- final_data %>%
  select(-one_of(columnas_a_eliminar))

data_pii<-(DynamicCancerDriverKM::PPI)

data_piin <- data_pii %>%
  pivot_longer(cols = c(`Input-node Gene Symbol`, `Output-node Gene Symbol`), names_to = "variable", values_to = "gen") %>%
  group_by(gen, variable) %>%
  summarise(frecuencia = n()) %>%
  pivot_wider(names_from = variable, values_from = frecuencia, values_fill = 0)

data_piinR <- data_piin %>%
  mutate(total_mode = `Input-node Gene Symbol` + `Output-node Gene Symbol`) %>%
  select(total_mode) %>%
  arrange(desc(total_mode))

# Mostrar el resultado
print(data_piinR)

genes_en_final_data <- colnames(final_data_filtrado)

# Filtrar data_piinR para mantener solo los genes que están en final_data_filtrado
data_piinR_filtrado <- data_piinR %>%
  filter(gen %in% genes_en_final_data)

# Mostrar el resultado
print(data_piinR_filtrado)

view(genes_en_final_data)
