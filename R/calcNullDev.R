calcNullDev <- function(X, y, group, family) {
  form <- if (any(group==0)) formula(y~X[, group==0]) else formula(y~1)
  fit <- glm(form, family=family)
  mean(loss.grpreg(y, predict(fit, type="response"), family))
}
