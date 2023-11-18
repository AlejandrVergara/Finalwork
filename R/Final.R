library(tidyverse)
library(caret)
library(class)
library(gmodels)
library(psych)
library(DynamicCancerDriverKM)

datanormal <-(DynamicCancerDriverKM::BRCA_normal)
dataPt <-(DynamicCancerDriverKM::BRCA_PT)
final_data <- bind_rows(datanormal, dataPt)


porcentaje_menor_10 <- final_data %>%
  summarise_all(~ mean(. <600, na.rm = TRUE))


columnas_a_eliminar <- names(porcentaje_menor_10[, porcentaje_menor_10 >= 0.8])


final_data_filtrado <- final_data %>%
  select(-one_of(columnas_a_eliminar))

final_data_filtrado2 <- final_data_filtrado

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


print(data_piinR)

final_data_filtradox<-colnames(final_data_filtrado)[ 8:ncol(final_data_filtrado)]
aux2 <- AMCBGeneUtils::changeGeneId(final_data_filtradox, from = "Ensembl.ID")

names(final_data_filtrado)[8:ncol(final_data_filtrado)] <- aux2$HGNC.symbol


genes_en_final_data <- colnames(final_data_filtrado)


data_piinR_filtrado <- data_piinR %>%
  filter(gen %in% genes_en_final_data)

Predictores <- as.vector(head(data_piinR_filtrado[, 1], 100))

Predictores <- as.character(unlist(Predictores))

colnames(final_data_filtrado)[is.na(colnames(final_data_filtrado))] <- paste0("gen", seq_along(colnames(final_data_filtrado) == ""))
set.seed(1)
final_data_filtradoe <- final_data_filtrado %>%
  group_by(sample_type) %>%
  sample_n(123, replace = TRUE) %>%
  ungroup()

sample.index <- sample(1:nrow(final_data_filtradoe)
                       ,nrow(final_data_filtradoe)*0.7
                       ,replace = F)

train.data <- final_data_filtradoe[sample.index, c(Predictores, "sample_type"), drop = FALSE]
test.data <- final_data_filtradoe[-sample.index, c(Predictores, "sample_type"), drop = FALSE]


train.data$sample_type <- factor(train.data$sample_type)
test.data$sample_type <- factor(test.data$sample_type)

# Train the k-NN model
ctrl <- trainControl(method = "cv", p = 0.7)
knnFit <- train(sample_type ~ .
                , data = train.data
                , method = "knn", trControl = ctrl
                , preProcess = c("range") # c("center", "scale") for z-score
                , tuneLength = 50)

plot(knnFit)

# Make predictions
knnPredict <- predict(knnFit, newdata = test.data)

# Creates the confusion matrix
confusionMatrix(data = knnPredict, reference = test.data$sample_type)


#regresion lineal
final_data_filtradoe <- final_data_filtradoe %>%
  mutate(sample_type = ifelse(sample_type == "Solid Tissue Normal", 0, 1))

train.data <- final_data_filtradoe[sample.index, c(Predictores, "sample_type"), drop = FALSE]
test.data <- final_data_filtradoe[-sample.index, c(Predictores, "sample_type"), drop = FALSE]

ins_model <- lm(sample_type ~ ., data = train.data)

summary(ins_model)


# Train the model
train.control <- trainControl(method = "cv", number = 10 )
model <- train(sample_type ~ ., data = train.data, method = "lm",
               trControl = train.control)

# Summarize the results

print(model)
