library(lightgbm)
data(agaricus.train, package='lightgbm')
train <- agaricus.train
dtrain <- lgb.Dataset(train$data, label = train$label)
model <- lgb.train(
  params = list(
    objective = "regression"
    , metric = "l2"
  )
  , data = dtrain
)

data(agaricus.test, package='lightgbm')
test <- agaricus.test
dtest <- lgb.Dataset(test$data, label = test$label)
model
predictions <- predict(model,test$data)

data.frame(test$label, round(predictions))
InformationValue::confusionMatrix(test$label, predictions)