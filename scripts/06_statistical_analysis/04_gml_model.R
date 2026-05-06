# ============================================================
# PACKAGES
# ============================================================
library(tidyverse)
library(readxl)
library(stringr)
library(broom)
library(scales)
library(lme4)
library(broom.mixed)

# ============================================================
# 1) LOAD FQ CONSUMPTION (EUROPE ONLY)
# ============================================================
fq <- read_excel("fq_long.xlsx") %>%
  rename(Country = Country,
         Year = Year,
         FQ_DDD = FQ_DDD) %>%
  filter(!is.na(Year), !is.na(FQ_DDD))

eu_countries <- fq %>% distinct(Country) %>% pull()

# ============================================================
# 2) LOAD CAMPYLOBACTER DATA AND CLEAN
# ============================================================
campy_raw <- read_excel("cleaned_campy_data_deduplicated.xlsx")

# Fix common issues automatically
campy <- campy_raw %>%
  # Fix the typo if present
  rename(Reistance = any_of("Reistance"),
         Reistance = any_of("Resistance")) %>%
  # Harmonize country names
  mutate(Country = recode(Country,
                          "Czech Republic" = "Czechia",
                          "UK" = "United Kingdom",
                          "England" = "United Kingdom",
                          "Scotland" = "United Kingdom",
                          "Wales" = "United Kingdom"
  )) %>%
  mutate(Year = suppressWarnings(as.integer(Year))) %>%
  filter(Country %in% eu_countries) %>%
  filter(!is.na(Country), !is.na(Year), !is.na(Reistance))

# ============================================================
# 3) COUNTRY–YEAR AGGREGATION
# ============================================================
campy_cy <- campy %>%
  group_by(Country, Year) %>%
  summarise(
    N_isolates = n(),
    Resistant  = sum(Reistance),
    Resistance_Prevalence = Resistant / N_isolates,
    .groups = "drop"
  ) %>%
  filter(N_isolates >= 10)

# ============================================================
# 4) MERGE WITH FQ CONSUMPTION + LAGS
# ============================================================
dat <- campy_cy %>%
  inner_join(fq, by = c("Country","Year")) %>%
  arrange(Country, Year) %>%
  group_by(Country) %>%
  mutate(
    FQ_DDD_lag1 = lag(FQ_DDD, 1)
  ) %>%
  ungroup() %>%
  mutate(
    FQ_scaled      = FQ_DDD,       # OR per +1 DDD
    FQ_lag1_scaled = FQ_DDD_lag1
  )

# ============================================================
# 5) OUTBREAK-SAMPLING CORRECTION
# ============================================================
cap_threshold <- quantile(dat$N_isolates, 0.95, na.rm = TRUE)

dat <- dat %>%
  mutate(
    N_capped = pmin(N_isolates, cap_threshold),
    log_weight = log1p(N_isolates)
  )

# ============================================================
# 6) MIXED-EFFECTS MODELS
# ============================================================

# ---- Current-year model ----
glmer_curr <- glmer(
  cbind(Resistant, N_isolates - Resistant) ~ FQ_scaled + (1 | Country),
  data = dat,
  family = binomial(),
  weights = N_capped
)

# ---- Lag-1 model ----
dat_lag1 <- dat %>% filter(!is.na(FQ_lag1_scaled))

glmer_lag1 <- glmer(
  cbind(Resistant, N_isolates - Resistant) ~ FQ_lag1_scaled + (1 | Country),
  data = dat_lag1,
  family = binomial(),
  weights = N_capped
)

# ============================================================
# 7) EXTRACT ORs PER +1 DDD
# ============================================================
or_from_glmer <- function(fit, term){
  smry <- tidy(fit, effects = "fixed")
  row  <- smry %>% filter(term == !!term)
  beta <- row$estimate
  se   <- row$std.error
  z    <- 1.96
  tibble(
    term = term,
    beta = beta,
    se = se,
    p = row$p.value,
    OR_per_1DDD = exp(beta),
    OR_low      = exp(beta - z*se),
    OR_high     = exp(beta + z*se)
  )
}

tab_curr <- or_from_glmer(glmer_curr, "FQ_scaled") %>%
  mutate(Model = "Mixed-effects (current-year)")

tab_lag1 <- or_from_glmer(glmer_lag1, "FQ_lag1_scaled") %>%
  mutate(Model = "Mixed-effects (lag-1)")

res_table <- bind_rows(tab_curr, tab_lag1)
print(res_table)

# ============================================================
# 8) DESCRIPTIVE COUNTRY-MEANS ANALYSIS
# ============================================================
country_means <- dat %>%
  group_by(Country) %>%
  summarise(
    FQ_DDD_mean = weighted.mean(FQ_scaled, w = N_isolates, na.rm = TRUE),
    ResPrev_mean = weighted.mean(Resistance_Prevalence, w = N_isolates, na.rm = TRUE),
    N_total = sum(N_isolates),
    n_years = n(),
    .groups = "drop"
  )

wls_fit <- lm(ResPrev_mean ~ FQ_DDD_mean, data = country_means, weights = N_total)
wls_tab <- tidy(wls_fit)
cor_test <- cor.test(country_means$FQ_DDD_mean, country_means$ResPrev_mean)

print(wls_tab)
print(cor_test)

























# ============================
# Packages
# ============================
library(tidyverse)
library(readxl)
library(stringr)
library(broom)
library(broom.mixed)
library(lme4)
library(performance)   # for check_overdispersion, R2, etc.
library(scales)

# ============================
# 1) Load data
# ============================

# Fluoroquinolone consumption (already long format: fq_long.xlsx)
fq <- read_excel("fq_long.xlsx") %>%
  janitor::clean_names() %>%
  rename(
    Country = country,
    Year    = year,
    FQ_DDD  = fq_ddd
  ) %>%
  filter(!is.na(Country), !is.na(Year), !is.na(FQ_DDD))

# Campylobacter dataset (deduplicated)
campy_raw <- read_excel("cleaned_campy_data_deduplicated.xlsx")

# Expecting columns like: Country, Year, Resistance (0/1) or Reistance
# Adjust this line if your column is named differently
campy <- campy_raw %>%
  janitor::clean_names() %>%
  mutate(
    year    = suppressWarnings(as.integer(year)),
    Country = country
  ) %>%
  # keep only rows with defined country, year, and resistance
  filter(!is.na(Country), !is.na(year), !is.na(resistance)) %>%
  rename(
    Year       = year,
    Resistant  = resistance
  )

# ============================
# 2) Aggregate to country–year
# ============================

campy_cy <- campy %>%
  group_by(Country, Year) %>%
  summarise(
    N_isolates  = n(),
    Resistant   = sum(Resistant),
    ResPrev     = Resistant / N_isolates,
    .groups = "drop"
  ) %>%
  # quality filter: avoid tiny cells that are unstable
  filter(N_isolates >= 10)

# Merge with FQ consumption
dat <- campy_cy %>%
  inner_join(fq, by = c("Country", "Year")) %>%
  arrange(Country, Year) %>%
  group_by(Country) %>%
  mutate(
    FQ_DDD_lag1 = dplyr::lag(FQ_DDD, 1)
  ) %>%
  ungroup() %>%
  # scale predictors to avoid insane coefficients
  mutate(
    FQ_scaled      = as.numeric(scale(FQ_DDD)),
    FQ_lag1_scaled = as.numeric(scale(FQ_DDD_lag1))
  )

# Keep a lag-1 subset
dat_lag1 <- dat %>% filter(!is.na(FQ_lag1_scaled))

# ============================
# 3) Helper: model summary table
# ============================

summarise_model <- function(fit, name, term){
  # term: predictor name to extract (e.g. "FQ_scaled" or "FQ_lag1_scaled")
  # Works for glm and glmer
  if (inherits(fit, "glmerMod")) {
    tt <- broom.mixed::tidy(fit, effects = "fixed")
    gl <- broom.mixed::glance(fit)
  } else {
    tt <- broom::tidy(fit)
    gl <- broom::glance(fit)
  }
  
  row <- tt %>% filter(term == !!term)
  if (nrow(row) == 0) {
    return(tibble(
      Model = name,
      term  = term,
      beta  = NA_real_,
      se    = NA_real_,
      p     = NA_real_,
      OR    = NA_real_,
      OR_low = NA_real_,
      OR_high = NA_real_,
      AIC   = gl$AIC,
      BIC   = gl$BIC,
      overdispersion = NA_real_,
      R2_marginal = NA_real_,
      R2_conditional = NA_real_
    ))
  }
  
  beta <- row$estimate
  se   <- row$std.error
  z    <- 1.96
  
  # OR per 1 SD of FQ (since we scaled)
  OR    <- exp(beta)
  OR_lo <- exp(beta - z*se)
  OR_hi <- exp(beta + z*se)
  
  # Overdispersion and R2 (for glmer only; for glm we approximate)
  if (inherits(fit, "glmerMod")) {
    od  <- tryCatch(performance::check_overdispersion(fit)$ratio, error = function(e) NA_real_)
    r2  <- tryCatch(performance::r2_nakagawa(fit), error = function(e) NULL)
    R2m <- if (!is.null(r2)) r2$R2_marginal else NA_real_
    R2c <- if (!is.null(r2)) r2$R2_conditional else NA_real_
  } else {
    # crude overdispersion for glm
    rdf <- df.residual(fit)
    od  <- sum(residuals(fit, type = "pearson")^2) / rdf
    R2m <- NA_real_
    R2c <- NA_real_
  }
  
  tibble(
    Model = name,
    term  = term,
    beta  = beta,
    se    = se,
    p     = row$p.value,
    OR    = OR,
    OR_low = OR_lo,
    OR_high = OR_hi,
    AIC   = gl$AIC,
    BIC   = gl$BIC,
    overdispersion = od,
    R2_marginal = R2m,
    R2_conditional = R2c
  )
}

# ============================
# 4) Fit competing models
# ============================

# 4.1 Simple weighted GLM (current-year only)
m_glm_curr <- glm(
  cbind(Resistant, N_isolates - Resistant) ~ FQ_scaled,
  data = dat,
  family = binomial(),
  weights = N_isolates
)

# 4.2 Mixed-effects: random intercept for Country (current-year)
m_mixed_country <- glmer(
  cbind(Resistant, N_isolates - Resistant) ~ FQ_scaled + (1 | Country),
  data = dat,
  family = binomial(),
  weights = N_isolates,
  control = glmerControl(optimizer = "bobyqa")
)

# 4.3 Mixed-effects: random intercepts for Country and Year
m_mixed_country_year <- glmer(
  cbind(Resistant, N_isolates - Resistant) ~ FQ_scaled + (1 | Country) + (1 | Year),
  data = dat,
  family = binomial(),
  weights = N_isolates,
  control = glmerControl(optimizer = "bobyqa")
)

# 4.4 Mixed-effects with lag-1 (Country + Year)
m_mixed_lag1 <- glmer(
  cbind(Resistant, N_isolates - Resistant) ~ FQ_scaled + FQ_lag1_scaled +
    (1 | Country) + (1 | Year),
  data = dat_lag1,
  family = binomial(),
  weights = N_isolates,
  control = glmerControl(optimizer = "bobyqa")
)

# ============================
# 5) Compare models
# ============================

res_models <- bind_rows(
  summarise_model(m_glm_curr,           "GLM weighted (current-year)", "FQ_scaled"),
  summarise_model(m_mixed_country,      "Mixed-effects: Country",      "FQ_scaled"),
  summarise_model(m_mixed_country_year, "Mixed-effects: Country+Year", "FQ_scaled"),
  summarise_model(m_mixed_lag1,         "Mixed-effects: lag-1 (FQ_lag1_scaled)", "FQ_lag1_scaled")
) %>%
  mutate(
    beta      = round(beta, 3),
    se        = round(se, 3),
    OR        = round(OR, 3),
    OR_low    = round(OR_low, 3),
    OR_high   = round(OR_high, 3),
    AIC       = round(AIC, 1),
    BIC       = round(BIC, 1),
    overdispersion = round(overdispersion, 2),
    R2_marginal    = round(R2_marginal, 3),
    R2_conditional = round(R2_conditional, 3),
    p         = signif(p, 3)
  )

print(res_models)

# ============================
# 6) Optional: country-mean WLS + correlation
# ============================

country_means <- dat %>%
  group_by(Country) %>%
  summarise(
    FQ_DDD_mean = weighted.mean(FQ_DDD, w = N_isolates, na.rm = TRUE),
    ResPrev_mean = weighted.mean(ResPrev, w = N_isolates, na.rm = TRUE),
    N_total = sum(N_isolates),
    n_years = n(),
    .groups = "drop"
  )

wls_fit <- lm(ResPrev_mean ~ FQ_DDD_mean, data = country_means, weights = N_total)
wls_tab <- tidy(wls_fit)
cor_test <- cor.test(country_means$FQ_DDD_mean, country_means$ResPrev_mean, method = "pearson")

print(wls_tab)
print(cor_test)

# ============================
# 7) Optional: quick plot
# ============================

ggplot(country_means,
       aes(x = FQ_DDD_mean, y = ResPrev_mean, size = N_total)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, aes(weight = N_total)) +
  scale_size_continuous(range = c(2, 10), name = "Total isolates") +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    x = "Average fluoroquinolone consumption (DDD per 1000 inhabitants/day)",
    y = "Average predicted FQ-resistant Campylobacter prevalence",
    title = "Europe: ecological association between FQ use and gyrA-based resistance"
  ) +
  theme_minimal(base_size = 12)

