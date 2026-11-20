# -----------------------------------------------------------------------------------------------------------------------------------------#

# Using restricted mean survival time and Cox regression to teach survival analysis and epidemiological bias: a Shakespearean worked example
# Dataset: Shakespeare_Dataset.xlsx
# Notes:
# # - The main denominator for incidence rates is character-acts, analogous to person-time in epidemiology.

# -----------------------------------------------------------------------------------------------------------------------------------------#

# Data, packages and parameters ----
data_file_name <- "Shakespeare_Dataset.xlsx" 
out_dir <- "your_path\\SHAKESPEARE" # put your own directory here
install_missing_packages <- TRUE # control whether to automatically install missing packages 
set.seed(20260601) # set seed for reproducibility
required_packages <- c("xlsx", "sysfonts", "performance", "showtext", "ggrepel", "readxl", "dplyr", "tidyr", "stringr", "survminer", "rnaturalearth", "rnaturalearthdata", "sf", "car",
                       "forcats", "janitor", "readr", "ggplot2", "broom", "sandwich", "lmtest", "survival", "MASS", "purrr", "survRM2", "patchwork")
if (install_missing_packages) { # if missing packages, here they installed from CRAN
  missing_packages <- required_packages[
    !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
  ]
  if (length(missing_packages) > 0) {
    install.packages(missing_packages, dependencies = TRUE)
  }
}
invisible(lapply(required_packages, library, character.only = TRUE))
# set some macro aspects of plots
#font_add_google("MedievalSharp", "medieval")
showtext_auto()
font2 <- "medieval" # set this font for subsequent plots in a more "Shakespearean" aesthetics
font2 <- "calibri"
bg_color <- "#fbf8f4"; par(bg = bg_color, family = font2)
showtext_auto()

# Paths and folders ----
# if not existing, create folders as stated in out_dir
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
# mark some tidy folders in which to put results
tables_dir <- file.path(out_dir, "R_outputs", "tables")
figures_dir <- file.path(out_dir, "R_outputs", "figures")
models_dir <- file.path(out_dir, "R_outputs", "models")
objects_dir <- file.path(out_dir, "R_outputs", "r_objects")
logs_dir <- file.path(out_dir, "R_outputs", "logs")

dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(models_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(objects_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

candidate_data_paths <- c( data_file_name, file.path(out_dir, data_file_name) )
data_path <- candidate_data_paths[file.exists(candidate_data_paths)][1]

if (is.na(data_path)) {
  stop(
    "Input file not found. Please place '", data_file_name,
    "' either in the current working directory or in:\n", out_dir
  )
}

message("Using dataset: ", normalizePath(data_path, mustWork = TRUE))
message("Saving outputs to: ", normalizePath(out_dir, mustWork = TRUE))

# Helper functions to load in Global Environment ----
# dichotomous var converter to 0,1
yn_to01 <- function(x) {
  x_chr <- str_squish(as.character(x))
  dplyr::case_when(
    is.na(x_chr) ~ NA_real_,
    str_to_lower(x_chr) %in% c("yes", "y", "true", "1", "si", "sì") ~ 1,
    str_to_lower(x_chr) %in% c("no", "n", "false", "0") ~ 0,
    TRUE ~ NA_real_
  )
}
# factor converter
safe_factor <- function(x, unknown_label = "Unknown") {
  x_chr <- str_squish(as.character(x))
  x_chr[is.na(x_chr) | x_chr == ""] <- unknown_label
  factor(x_chr)
}
# assign reference level to categorical variables
relevel_if_present <- function(x, ref) {
  x <- factor(x)
  if (ref %in% levels(x)) {
    stats::relevel(x, ref = ref)
  } else {
    x
  }
}
# calculate Poisson confidence intervals
poisson_ci_rate <- function(events, exposure, multiplier = 100) {
  if (is.na(events) || is.na(exposure) || exposure <= 0) {
    return(c(NA_real_, NA_real_))
  }
  stats::poisson.test(events, T = exposure)$conf.int * multiplier
}
# call poisson rate calculator above and add column to df
add_rate_columns <- function(df, event_col = "deaths", exposure_col = "character_acts_total") {
  df %>%
    rowwise() %>%
    mutate(
      cfr_percent = ifelse(n_characters > 0, 100 * .data[[event_col]] / n_characters, NA_real_),
      incidence_rate_per_100_character_acts =
        ifelse(.data[[exposure_col]] > 0, 100 * .data[[event_col]] / .data[[exposure_col]], NA_real_),
      incidence_rate_lower_95 =
        poisson_ci_rate(.data[[event_col]], .data[[exposure_col]], multiplier = 100)[1],
      incidence_rate_upper_95 =
        poisson_ci_rate(.data[[event_col]], .data[[exposure_col]], multiplier = 100)[2]
    ) %>%
    ungroup()
}
# robust Poisson model fit
fit_poisson_robust <- function(formula, data, cluster_var = "play_title") {
  response_var <- all.vars(formula)[1]
  model_vars <- unique(c(all.vars(formula), cluster_var))
  model_vars <- model_vars[model_vars %in% names(data)]
  
  dat <- data %>%
    dplyr::select(dplyr::all_of(model_vars)) %>%
    dplyr::filter(dplyr::if_all(dplyr::everything(), ~ !is.na(.)))
  
  if (nrow(dat) == 0) {
    stop("No complete cases available for model: ", deparse(formula))
  }
  
  model <- stats::glm(formula, family = poisson(link = "log"), data = dat)
  
  cluster <- dat[[cluster_var]]
  
  if (!is.null(cluster) && length(unique(cluster)) > 1) {
    vcov_robust <- sandwich::vcovCL(model, cluster = cluster, type = "HC0")
  } else {
    vcov_robust <- sandwich::vcovHC(model, type = "HC0")
  }
  
  ct <- lmtest::coeftest(model, vcov. = vcov_robust)
  
  results <- tibble::tibble(
    term = rownames(ct),
    estimate = as.numeric(ct[, 1]),
    robust_se = as.numeric(ct[, 2]),
    statistic = as.numeric(ct[, 3]),
    p_value = as.numeric(ct[, 4]),
    irr = exp(estimate),
    conf_low = exp(estimate - 1.96 * robust_se),
    conf_high = exp(estimate + 1.96 * robust_se),
    n_obs = nrow(dat),
    n_events = sum(dat[[response_var]], na.rm = TRUE)
  )
  
  dispersion <- sum(stats::residuals(model, type = "pearson")^2) / model$df.residual
  
  list(
    model = model,
    robust_table = results,
    model_data = dat,
    dispersion = dispersion,
    response_var = response_var
  )
}

# Import and cleaning ----
## import data ----
raw <- readxl::read_excel(data_path) %>% janitor::clean_names()
ncol(raw)
names(raw)
## make analysis_data as an operative ds ----
analysis_data <- raw %>%
  mutate(
    # variable format
    character_id = as.integer(character_id),
    character_name = str_squish(as.character(character_name)),
    play_title = str_squish(as.character(play_title)),
    genre = str_squish(as.character(genre)),
    profession_primary = str_squish(as.character(profession_primary)),
    social_status = str_squish(as.character(social_status)),
    gender = str_squish(as.character(gender)),
    character_status = str_squish(as.character(character_status)),
    acts_survived = as.numeric(acts_survived),
    death_type_clean = str_squish(as.character(death_type)),
    death_on_stage_clean = str_squish(as.character(death_on_stage)),
    supernatural_involvement_clean = str_squish(as.character(supernatural_involvement)),
    betrayal_exposure_clean = str_squish(as.character(betrayal_exposure)),
    # dichotomise
    death_event = yn_to01(dies_during_play),
    suicide_event = yn_to01(suicide),
    poisoning_event = yn_to01(poisoning),
    execution_event = yn_to01(execution),
    battle_related_event = yn_to01(battle_related),
    messenger_role_event = yn_to01(messenger_role),
    supernatural_exposure = ifelse(
      is.na(supernatural_involvement_clean) | supernatural_involvement_clean == "" | str_to_lower(supernatural_involvement_clean) == "none", 0, 1)
  )
# replace NAs with 0
analysis_data <- analysis_data %>% mutate(
  death_event = ifelse(is.na(death_event), 0, death_event),
  suicide_event = ifelse(is.na(suicide_event), 0, suicide_event),
  poisoning_event = ifelse(is.na(poisoning_event), 0, poisoning_event),
  execution_event = ifelse(is.na(execution_event), 0, execution_event),
  battle_related_event = ifelse(is.na(battle_related_event), 0, battle_related_event),
  messenger_role_event = ifelse(is.na(messenger_role_event), 0, messenger_role_event)
)
analysis_data <- analysis_data %>% mutate(
  # Character-acts are the denominator for incidence analyses.
  # Values <=0 are set to 0.5 to avoid zero exposure in offset models.
  # This mainly affects prologues/choruses/liminal figures, if present.
  character_acts = dplyr::case_when(
    is.na(acts_survived) ~ NA_real_,
      acts_survived <= 0 ~ 0.5,
      TRUE ~ acts_survived
    ),
    year_cut = cut(year, breaks = c(1500, 1600, 1608, Inf), labels = c("1588-1600", "1601-1608", "1609-1613"), right = T), # categorised in Early period 1600, Middle period 1601-1608, Late period 1609-1623
    death_on_stage_category = dplyr::case_when(
      death_event == 0 ~ "No death",
      str_detect(death_on_stage_clean, regex("^yes$", ignore_case = TRUE)) ~ "On-stage",
      str_detect(death_on_stage_clean, regex("off", ignore_case = TRUE)) ~ "Off-stage",
      str_detect(death_on_stage_clean, regex("^no$", ignore_case = TRUE)) ~ "Off-stage or not shown",
      TRUE ~ "Unknown"
    ),
    onstage_death_event = ifelse(death_on_stage_category == "On-stage", 1, 0),
    offstage_death_event = ifelse( death_event == 1 & death_on_stage_category %in% c("Off-stage", "Off-stage or not shown"), 1, 0 ),
    violent_death_event = ifelse( 
      death_event == 1 & str_detect( death_type_clean, regex("violent|battle|execution|murder|killed|stab|poison|suicide", ignore_case = TRUE)), 1, 0 ),
    natural_death_event = ifelse(
      death_event == 1 & str_detect(death_type_clean, regex("natural", ignore_case = TRUE)), 1, 0 ),
    supernatural_death_event = ifelse(death_event == 1 & str_detect(death_type_clean, regex("supernatural", ignore_case = TRUE)), 1, 0 ),
    betrayal_exposure_cat = dplyr::case_when( str_to_lower(betrayal_exposure_clean) %in% c("low", "medium", "high") ~ betrayal_exposure_clean, TRUE ~ "Unknown")
)

# Profession classification and variable merging for modelling purposes
analysis_data <- analysis_data %>%
  mutate(
    profession_lower = str_to_lower(profession_primary),
    profession_group = dplyr::case_when(
      str_detect(profession_lower, "king|queen|prince|princess|duke|duchess|emperor|monarch|ruler|governor") ~
        "Royalty and rulers",
      str_detect(profession_lower, "noble|lord|lady|gentleman|gentlewoman|courtier|count|earl|baron|knight|senator|tribune|counsellor|councillor|advisor|ambassador") ~
        "Nobility and court",
      str_detect(profession_lower, "soldier|general|captain|officer|guard|warrior|marshal|military") ~
        "Military and security",
      str_detect(profession_lower, "messenger|herald|chorus|prologue|epilogue") ~
        "Messengers",
      str_detect(profession_lower, "friar|priest|clergy|cardinal|bishop|nun|abbess|prophet|soothsayer|curate|parish") ~
        "Clergy and prophecy",
      str_detect(profession_lower, "physician|doctor|healer|apothecary|surgeon") ~
        "Medical professions",
      str_detect(profession_lower, "servant|attendant|page|steward|nurse|waiting|butler|housekeeper|jailer|gaoler|executioner") ~
        "Servants and institutional staff",
      str_detect(profession_lower, "witch|ghost|spirit|god|goddess|apparition|fairy|supernatural") ~
        "Supernatural beings",
      str_detect(profession_lower, "merchant|tradesman|artisan|fisherman|sailor|mariner|pirate|shepherd|clown|fool|jester|citizen|commoner|laborer|rogue|bawd|courtesan|host|innkeeper|musician|singer") ~
        "Commoners, trades, and performers",
      TRUE ~ "Other or unclassified"
    ),
    genre = safe_factor(genre),
    gender = safe_factor(gender),
    social_status = safe_factor(social_status),
    character_status = safe_factor(character_status),
    profession_group = factor(profession_group),
    supernatural_exposure = factor(supernatural_exposure),
    # need to merge betrayal exposure categories low and medium
    betrayal_exposure_cat = ifelse(betrayal_exposure_cat %in% "High", "High", "Low/Medium"),
    betrayal_exposure_cat = factor(betrayal_exposure_cat, levels = c("Low/Medium", "High")),
    # need to merge character status
    character_status = ifelse(character_status %in% "Minor", "Minor", ifelse(character_status %in% c("Primary antagonist", "Primary protagonist", "Supporting antagonist"), "Protagonist or antagonist", "Supporting")),
    messenger_role_factor = factor( ifelse(messenger_role_event == 1, "Messenger role", "No messenger role") )
  ) %>%
  mutate(
    profession_group_model = forcats::fct_lump_min(profession_group, min = 10, other_level = "Other or rare professions"),
    profession_group_model = relevel_if_present(profession_group_model, "Commoners, trades, and performers"),
    gender = relevel_if_present(gender, "Man"),
    social_status = relevel_if_present(social_status, "Commoner"),
    character_status = relevel_if_present(character_status, "Minor"),
    supernatural_exposure = relevel_if_present(supernatural_exposure, "None"),
    messenger_role_factor = relevel_if_present(messenger_role_factor, "No messenger role")
  )


# Descriptive analyses ----
# extract some summary stats overall
overall_summary <- analysis_data %>%
  summarise(
    n_characters = n(),
    n_plays = n_distinct(play_title),
    n_genres = n_distinct(genre),
    character_acts_total = sum(character_acts, na.rm = TRUE),
    deaths = sum(death_event, na.rm = TRUE),
    onstage_deaths = sum(onstage_death_event, na.rm = TRUE),
    offstage_deaths = sum(offstage_death_event, na.rm = TRUE),
    violent_deaths = sum(violent_death_event, na.rm = TRUE),
    natural_deaths = sum(natural_death_event, na.rm = TRUE),
    supernatural_deaths = sum(supernatural_death_event, na.rm = TRUE),
    suicides = sum(suicide_event, na.rm = TRUE),
    poisonings = sum(poisoning_event, na.rm = TRUE),
    executions = sum(execution_event, na.rm = TRUE),
    battle_related = sum(battle_related_event, na.rm = TRUE)
  ) %>%
  add_rate_columns()
# extract some summary stats by genre, play, social status, and so on...
by_genre <- analysis_data %>%
  group_by(genre) %>%
  summarise(
    n_characters = n(),
    n_plays = n_distinct(play_title),
    character_acts_total = sum(character_acts, na.rm = TRUE),
    deaths = sum(death_event, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  add_rate_columns() %>%
  arrange(desc(incidence_rate_per_100_character_acts))

by_play <- analysis_data %>%
  group_by(play_title, genre) %>%
  summarise(
    n_characters = n(),
    character_acts_total = sum(character_acts, na.rm = TRUE),
    deaths = sum(death_event, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  add_rate_columns() %>%
  arrange(desc(incidence_rate_per_100_character_acts), play_title)

by_social <- analysis_data %>%
  group_by(social_status) %>%
  summarise(
    n_characters = n(),
    n_plays = n_distinct(play_title),
    character_acts_total = sum(character_acts, na.rm = TRUE),
    deaths = sum(death_event, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  add_rate_columns() %>%
  arrange(desc(incidence_rate_per_100_character_acts))

by_character <- analysis_data %>%
  group_by(character_status) %>%
  summarise(
    n_characters = n(),
    n_plays = n_distinct(play_title),
    character_acts_total = sum(character_acts, na.rm = TRUE),
    deaths = sum(death_event, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  add_rate_columns() %>%
  arrange(desc(incidence_rate_per_100_character_acts))

by_battle <- analysis_data %>%
  group_by(battle_related) %>%
  summarise(
    n_characters = n(),
    n_plays = n_distinct(play_title),
    character_acts_total = sum(character_acts, na.rm = TRUE),
    deaths = sum(death_event, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  add_rate_columns() %>%
  arrange(desc(incidence_rate_per_100_character_acts))

by_year <- analysis_data %>%
  group_by(year_cut) %>%
  summarise(
    n_characters = n(),
    n_plays = n_distinct(play_title),
    character_acts_total = sum(character_acts, na.rm = TRUE),
    deaths = sum(death_event, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  add_rate_columns() %>%
  arrange(desc(incidence_rate_per_100_character_acts))

by_betrayal <- analysis_data %>%
  group_by(betrayal_exposure_cat) %>%
  summarise(
    n_characters = n(),
    n_plays = n_distinct(play_title),
    character_acts_total = sum(character_acts, na.rm = TRUE),
    deaths = sum(death_event, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  add_rate_columns() %>%
  arrange(desc(incidence_rate_per_100_character_acts))

by_supernatural <- analysis_data %>%
  group_by(supernatural_exposure) %>%
  summarise(
    n_characters = n(),
    n_plays = n_distinct(play_title),
    character_acts_total = sum(character_acts, na.rm = TRUE),
    deaths = sum(death_event, na.rm = TRUE),
    onstage_deaths = sum(onstage_death_event, na.rm = TRUE),
    offstage_deaths = sum(offstage_death_event, na.rm = TRUE),
    violent_deaths = sum(violent_death_event, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  add_rate_columns() %>%
  arrange(desc(incidence_rate_per_100_character_acts))

by_profession_group <- analysis_data %>%
  group_by(profession_group) %>%
  summarise(
    n_characters = n(),
    n_plays = n_distinct(play_title),
    character_acts_total = sum(character_acts, na.rm = TRUE),
    deaths = sum(death_event, na.rm = TRUE),
    onstage_deaths = sum(onstage_death_event, na.rm = TRUE),
    offstage_deaths = sum(offstage_death_event, na.rm = TRUE),
    violent_deaths = sum(violent_death_event, na.rm = TRUE),
    suicides = sum(suicide_event, na.rm = TRUE),
    poisonings = sum(poisoning_event, na.rm = TRUE),
    executions = sum(execution_event, na.rm = TRUE),
    battle_related = sum(battle_related_event, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  add_rate_columns() %>%
  arrange(desc(incidence_rate_per_100_character_acts))

by_profession_primary <- analysis_data %>%
  group_by(profession_primary) %>%
  summarise(
    n_characters = n(),
    n_plays = n_distinct(play_title),
    character_acts_total = sum(character_acts, na.rm = TRUE),
    deaths = sum(death_event, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  add_rate_columns() %>%
  arrange(desc(deaths), desc(incidence_rate_per_100_character_acts), profession_primary)

by_country <- analysis_data %>%
  group_by(country) %>%
  summarise(
    n_characters = n(),
    n_plays = n_distinct(play_title),
    character_acts_total = sum(character_acts, na.rm = TRUE),
    deaths = sum(death_event, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  add_rate_columns() %>%
  arrange(desc(deaths), desc(incidence_rate_per_100_character_acts), country)

# summary of death counts
death_modes <- analysis_data %>%
  filter(death_event == 1) %>%
  count(death_type_clean, death_on_stage_clean, sort = TRUE, name = "n_deaths")
death_modes <- death_modes %>% pivot_wider(names_from = death_on_stage_clean, values_from = n_deaths, values_fill = 0)
death_type_by_profession <- analysis_data %>%
  filter(death_event == 1) %>%
  count(profession_group, death_type_clean, sort = TRUE, name = "n_deaths") %>%
  group_by(profession_group) %>%
  mutate(percent_within_profession_deaths = 100 * n_deaths / sum(n_deaths)) %>%
  ungroup()

# % of all deaths
table(analysis_data$death_type[analysis_data$death_event==1]) %>% prop.table()*100
# deaths per 100 character-acts
table(analysis_data$death_type[analysis_data$death_event==1])/sum(analysis_data$character_acts)*100
# write a xlsx file with death descriptives
write.xlsx(death_modes %>% as.data.frame(), file = file.path(tables_dir, "deaths-descriptive.xlsx"), sheetName = "Death-modes", row.names = F)
write.xlsx(death_type_by_profession %>% as.data.frame(), file = file.path(tables_dir, "deaths-descriptive.xlsx"), sheetName = "Death-type-by-profession", row.names = F, append = T)

# Modelling ----
# 1) Poisson models assume constant hazards across acts
# 2) Cox PH in a pseudo-discrete time setting such as this implies non-proportional hazards
# 3) Restricted Mean Survival Time
# 4) pooled logistic regression

## Poisson models: explore models by assuming constant hazard across acts (very unlikely) ----
# Outcome: death during the play.
# Offset: log(character_acts), analogous to person-time. To account for time at risk during plau
# Robust SEs: clustered by play_title.

# builf model formula
COV <- "social_status + genre + gender + character_status + supernatural_exposure + betrayal_exposure_cat + year_cut"
fom <- as.formula(paste0("death_event ~ ", COV, " + offset(log(character_acts))"))

# build models by adding variables
main_model_1 <- fit_poisson_robust( death_event ~ social_status + offset(log(character_acts)), data = analysis_data )
summary(main_model_1$model)
main_model_2 <- fit_poisson_robust( death_event ~ social_status + genre + gender + offset(log(character_acts)), data = analysis_data )
summary(main_model_2$model)
lrtest(main_model_1$model, main_model_2$model) # compare with model 1
main_model_3 <- fit_poisson_robust(fom, data = analysis_data)
summary(main_model_3$model)
car::vif(main_model_3$model) # VIF < 5 for genre and year_cut
lrtest(main_model_2$model, main_model_3$model) # compare with model 2

main_model_summary <- bind_rows(
  main_model_1$robust_table %>% mutate(model = "Model 1: profession only"),
  main_model_2$robust_table %>% mutate(model = "Model 2: + genre, gender, social status"),
  main_model_3$robust_table %>% mutate(model = "Model 3: fully adjusted")
) %>% dplyr::select(model, everything())

main_model_summary <- main_model_summary %>% mutate(
  IRR = sprintf("%.2f [%.2f, %.2f]", irr, conf_low, conf_high)
)

write.xlsx(main_model_summary, file.path(models_dir, "models.xlsx"))

# check overdispersion and negative binomial sensitivity 
overdispersion_summary <- tibble::tibble(
  model = c("main_model_1", "main_model_2", "main_model_3"),
  dispersion = c(main_model_1$dispersion, main_model_2$dispersion, main_model_3$dispersion)
)
overdispersion_summary

# Negative binomial sensitivity model.
# This is mainly a robustness check; convergence may be unstable with sparse categories.
nb_data <- analysis_data %>%
  dplyr::select( death_event, profession_group_model, genre, gender, social_status, character_status, messenger_role_factor, 
                 supernatural_exposure, betrayal_exposure_cat, character_acts, play_title, year_cut) %>%
  filter(if_all(everything(), ~ !is.na(.)))

nb_model <- MASS::glm.nb(fom, data = nb_data)
broom::tidy(nb_model, exponentiate = TRUE, conf.int = TRUE)
summary(nb_model)
AIC(nb_model, main_model_3$model) # compare with Poisson model 3

## Survival analysis ----
# be sure that character acts and death event is non-missing
survival_data <- analysis_data %>% filter(!is.na(character_acts), !is.na(death_event)) %>%
  mutate( surv_time = character_acts, # define survival time in character-acts
          surv_event = death_event, # define survival event = death (1) or censoring (0) at some point up to play end (time = 5)
          profession_group_surv = profession_group_model ) # just rename profession to be used later
### add the example of a FULLY SYNTHETIC (not based on actual Shakespearean plays), time-varying supernatural exposure ----
# randomly assign time to exposure among the exposed before death or censoring, assuming association with death
set.seed(2604) # set seed for reproducibility of CUSTOM-DEFINED TIME_VARYING SUPERNATURAL EXPOSURE
survival_data <- survival_data %>% mutate(
  sup_bin = supernatural_exposure == 1, # instrumental flag
  # start-time draws: skewed early for those who die, skewed late for censored
  t_early = sample(1:5, n(), replace = TRUE, prob = c(.40, .25, .20, .10, .05)), # randomly sample time 1 to 5 based on those arbitrary prior probabilities
  t_late  = sample(1:5, n(), replace = TRUE, prob = c(.05, .20, .25, .20, .30)), # different distribution to allow for assumed association betweeen death and time to supernatural exposure
  t_sup   = ifelse(sup_bin, ifelse(surv_event == 1, t_early, t_late), NA_real_),
  # exposure must start strictly before end of follow-up
  timeto_supernatural_exposure = ifelse(!is.na(t_sup) & t_sup >= surv_time, NA_real_, t_sup),
  # Inf = never exposed (for tt() approach)
  tv_supernatural_exposure = ifelse(is.na(timeto_supernatural_exposure), Inf, timeto_supernatural_exposure)
) %>% dplyr::select(-t_early, -t_late)
table(survival_data$tv_supernatural_exposure, survival_data$surv_time, exclude = F)


### draw some KM curves ----### draw some KM curves ----surv_time
km_by_genre <- survival::survfit( survival::Surv(surv_time, surv_event) ~ genre, data = survival_data )
km_by_profession <- survival::survfit( survival::Surv(surv_time, surv_event) ~ social_status, data = survival_data )
km_by_year <- survival::survfit( survival::Surv(surv_time, surv_event) ~ year_cut, data = survival_data )

col_genre <- c("History" = "#CCFBFE", "Romance" = "#363636", "Tragedy" = "#247BA0")
col_prof <- colorRampPalette( c("#0F8B8D", "#CCFBFE", "#143642", "#363636"))(length(km_by_profession$strata))

# by genre

png(file = paste0(figures_dir, "/KM.png"), width = 2000, height = 1000, res = 150)

par(mfrow = c(1, 3), bg = bg_color, family = font2, cex = 1.1, cex.lab = 1.2, cex.main = 1.2, cex.axis = 1.2)
plot(km_by_genre, ylim = c(.5, 1), lty = c(1, 2, 3), lwd = 2, col = col_genre, xlab = "Time (acts)", ylab = "Survival (%)", 
     mark.time = TRUE, mark = 3, cex = 0.6, bty = "l", axes = FALSE, main = "Survival by genre")
axis(1, at = seq(0, 5, by = 1))
axis(2, at = seq(0, 1, by = .2), las = 2)
abline(h = seq(0, 5, by = 1), col = "gray90", lty = 3)
abline(v = seq(0, 5, by = 1), col = "gray90", lty = 3)
legend("bottomleft", legend = str_replace(names(km_by_genre$strata),"genre=", ""), lty = c(1, 2, 3),
       lwd = 2, col = col_genre, bty = "n", cex = 0.7)

# by profession 
plot(km_by_profession, ylim = c(.5, 1), col = col_prof, lty = 2, lwd = 2, xlab = "Time (acts)", ylab = "Survival (%)", 
     mark.time = TRUE, mark = 3, cex = 0.6, bty = "l", axes = FALSE, main = "Survival by social status")
axis(1, at = seq(0, 5, by = 1))
axis(2, at = seq(0, 1, by = .2), las = 2)
abline(h = seq(0, 5, by = 1), col = "gray90", lty = 3)
abline(v = seq(0, 5, by = 1), col = "gray90", lty = 3)
legend("bottomleft", legend = str_replace(names(km_by_profession$strata),"social_status=", ""), lty = 2,
       lwd = 2, col = col_prof, bty = "n", cex = 0.7)

# by year 
plot(km_by_year, ylim = c(.5, 1), col = col_prof, lty = 2, lwd = 2, xlab = "Time (acts)", ylab = "Survival (%)", 
     mark.time = TRUE, mark = 3, cex = 0.6, bty = "l", axes = FALSE, main = "Survival by play calendar year")
axis(1, at = seq(0, 5, by = 1))
axis(2, at = seq(0, 1, by = .2), las = 2)
abline(h = seq(0, 5, by = 1), col = "gray90", lty = 3)
abline(v = seq(0, 5, by = 1), col = "gray90", lty = 3)
legend("bottomleft", legend = c("Early period (1588-1600)", "Middle period (1601-1608)", "Late period (1609-1613)"), lty = 2,
       lwd = 2, col = col_prof, bty = "n", cex = 0.7)
par(mfrow = c(1, 1))
dev.off()

## Cox models ----
taus <- c(3, 5) # at what time horizons I want to calculate Risk Differences
tt_fun <- function(x, t, ...) as.numeric(t > x) # function to handle time-varying covariate
ties = "breslow"
n_boot = 2000

### helper functions ----
# calculate risk difference (RD) from Cox models
cox_rd <- function(fit, data, group_var, taus, basefit) {
  basefit <- survival::survfit(fit)
  if (is.null(basefit)) return(hr_out) # HR only for tt() models (or if baseline curve fails for any other reason)
  #computed ONCE per bootstrap iteration (baseline survival, centered at covariate means)
  # here we only compute individual linear predictors and combine with the shared baseline -> no repeated survfit() calls
  levs <- levels(as.factor(data[[group_var]]))
  ref <- levs[1]
  time_idx <- purrr::map_int(taus, function(tau) max(which(basefit$time <= tau))) # index of each tau on the shared baseline time grid
  
  risk_at_taus <- function(lvl) {
    nd <- data
    nd[[group_var]] <- factor(lvl, levels = levs) # set everyone to this level, keep other covariates as observed
    lp <- predict(fit, newdata = nd, type = "lp") # centered linear predictor per subject (same centering as basefit)
    s0 <- basefit$surv[time_idx] # baseline survival at the two taus (shared across all group_vars/levels)
    # S_i(t) = S0(t)^exp(lp_i)  -> average over subjects, then 1 - S for cumulative incidence
    purrr::map_dbl(seq_along(taus), function(k) {
      mean(1 - s0[k]^exp(lp))
    })
  }
  
  purrr::map_dfr(levs, function(lvl) {
    tibble::tibble(level = lvl, tau = taus, risk = risk_at_taus(lvl))
  }) %>%
    dplyr::group_by(tau) %>%
    dplyr::mutate(risk_difference = risk - risk[level == ref]) %>%
    dplyr::ungroup() %>% mutate(term = paste0(group_var, level))
}
# Bootstrap helper for Cox models (cluster resample by play_title)
cox_boot_ci <- function(formula, data, tt_fun = NULL, ties = "breslow", n_boot = 2000, label = "") {
  play_ids <- unique(data$play_title) # take unique play titles
  taus <- c(3, 5) # at what time horizons I want to calculate Risk Differences
  rhs_terms <- attr(stats::terms(formula), "term.labels") # for what variables in the formula
  groupvars <- rhs_terms[!grepl("^cluster\\(|^tt\\(", rhs_terms)] # extract them
  
  # pre-split data by play_title once, outside the loop: turns 2000 x (filter per sampled play) into 2000 x (list lookup)
  data_by_play <- base::split(data, data$play_title) # list of plays' characters
  
  set.seed(2604) # reproducible, same seed convention as before
  boots <- purrr::map_dfr(seq_len(n_boot), function(i) { # for each bootstrap sample i run the following function:
    sampled_plays <- sample(play_ids, length(play_ids), replace = TRUE) # resample plays with replacement
    d_boot <- dplyr::bind_rows(data_by_play[sampled_plays]) # single bind instead of n_plays separate filters
    
    formula_i <- formula # make formula in the current new d_boot environment to avoid issues
    environment(formula_i) <- environment()
    
    fit <- tryCatch(
      if (!is.null(tt_fun)) coxph(formula_i, data = d_boot, ties = ties, tt = tt_fun) # try fit the model
      else coxph(formula_i, data = d_boot, ties = ties),
      error = function(e) NULL
    )
    if (is.null(fit)) return(NULL) # skip iteration if fit failed
    
    hr_out <- broom::tidy(fit, exponentiate = TRUE) %>% dplyr::mutate(iteration = i) # store results into boots container
    
    # calculate risk difference at 3rd and 5th act
    rd_out <- purrr::map_dfr(groupvars, function(gv) { # RD for every covariate in the formula
      tryCatch(
        cox_rd(fit, d_boot, group_var = gv, taus = taus) %>%
          dplyr::rename(level = 1) %>% # cox_rd's first column is named after group_var; standardize to "level"
          dplyr::mutate(group_var = gv, iteration = i),
        error = function(e) NULL
      )
    })
    
    # resulting ds, minimise information and n cols
    if (nrow(rd_out) == 0) return(hr_out) # if RD failed entirely, still keep HR results for this iteration
    rd_out %>% dplyr::mutate(term = paste0(group_var, level)) %>%
      dplyr::left_join(hr_out %>% dplyr::select(-iteration) %>% dplyr::rename(HR = estimate), by = "term")
  }, .progress = paste("Bootstrap:", label)) # progress bar
  boots # return boots container full
}
# add bootstrap CI (HR_boot) alongside the existing analytic-CI HR column
add_boot_hr <- function(res, tidy_full, boots) {
  quant <- function(x, p) stats::quantile(x, p, na.rm = TRUE, names = FALSE) # take quantiles despite NAs if any
  boot_ci <- boots %>% 
    dplyr::group_by(term) %>% # for each model term
    dplyr::summarise(HR_boot_lo = quant(HR, .025), HR_boot_hi = quant(HR, .975), 
                     RD_at_act3_boot_lo = quant(RD_at_act3, .025), RD_at_act3_boot_hi = quant(RD_at_act3, .975), 
                     RD_at_act5_boot_lo = quant(RD_at_act5, .025), RD_at_act5_boot_hi = quant(RD_at_act5, .975), 
                     .groups = "drop") # take only 0.25 and .975 quantiles
  tidy_full %>%  # order and format
    dplyr::select(term, HR, RD_at_act3, RD_at_act5) %>%
    dplyr::left_join(boot_ci, by = "term") %>%
    dplyr::mutate(HR_boot = sprintf("%.2f [%.2f, %.2f]", HR, HR_boot_lo, HR_boot_hi),
                  RD_at_act3_boot = sprintf("%.2f [%.2f, %.2f]", RD_at_act3, RD_at_act3_boot_lo, RD_at_act3_boot_hi),
                  RD_at_act5_boot = sprintf("%.2f [%.2f, %.2f]", RD_at_act5, RD_at_act5_boot_lo, RD_at_act5_boot_hi)) %>%
    dplyr::select(term, HR_boot, RD_at_act3_boot, RD_at_act5_boot) %>%
    dplyr::right_join(res, by = "term") %>%
    dplyr::relocate(term, HR_boot, RD_at_act3_boot, RD_at_act5_boot) %>% dplyr::select(-HR)
}

### adjusted models alternatives ----
# please, refer to the hypothesised DAG depicted in Figure A1 in Supplementary Materials
CRUDE <- "social_status"
COV_ENTRY <- "social_status + genre + gender + year_cut" # entry adjustment
ENTRY <- strsplit(COV_ENTRY, " \\+ ")[[1]]

COV_EXTENDED <- paste(
  COV_ENTRY,
  "character_status + supernatural_exposure + betrayal_exposure_cat",
  sep = " + "
) #extended descriptive model adjustment
EXTENDED <- strsplit(COV_EXTENDED, " \\+ ")[[1]]

#### cox crude ----
cox_crude <- coxph(
  Surv(surv_time, surv_event) ~ social_status + cluster(play_title),
  data = survival_data,
  ties = "breslow"
)
survival::cox.zph(cox_crude) # check PH assumption by formal test: whether scaled Schoenfeld residuals / effects of covariates are constant over time at risk
cox.crude.tidy <- broom::tidy(cox_crude, exponentiate = TRUE, conf.int = TRUE) %>% mutate(
  HR = sprintf("%.2f [%.2f, %.2f]", estimate, conf.low, conf.high)
) 
# add risk difference and format model results
rd <- cox_rd(cox_crude, survival_data, group_var = "social_status", taus)
cox.crude.tidy <- cox.crude.tidy %>% left_join(rd, by = "term") %>% rename(fHR = HR, HR = estimate) %>% 
  dplyr::select(-risk) %>% pivot_wider(names_from = tau, names_prefix = "RD_at_act", values_from = risk_difference) # pivot taus as columns
cox.crude.res <- cox.crude.tidy %>% dplyr::select(term, HR)
# bootstrap call: may take a few minutes
boots.crude <- cox_boot_ci(formula = Surv(surv_time, surv_event) ~ social_status + cluster(play_title), data = survival_data, label = "Cox_crude")
boots.crude <- boots.crude %>% dplyr::select(-risk) %>% pivot_wider(names_from = tau, names_prefix = "RD_at_act", values_from = risk_difference) # pivot taus as columns
cox.crude.res <- add_boot_hr(cox.crude.res, cox.crude.tidy, boots.crude)
cox.crude.res
write.xlsx(cox.crude.res %>% as.data.frame(), file.path(models_dir, "models.xlsx"), sheetName = "Cox_crude", append = T, row.names = F, showNA = F)

#### cox entry adjustment ----
cox_entry_formula <- as.formula(paste0(
  "Surv(surv_time, surv_event) ~ ", # survival is a bi-dimensional quantity: time at which event may be 1 (death) or 0 (censoring)
  COV_ENTRY,
  " + cluster(play_title)"
))
cox_entry <- coxph(
  cox_entry_formula,
  data = survival_data, 
  ties = "breslow" # deal with ties in survival time with Breslow method
)
survival::cox.zph(cox_entry) # check PH assumption by formal test
ggcoxzph(cox.zph(cox_entry)) # Schoenfeld residuals plots
cox.entry.tidy <- broom::tidy(cox_entry, exponentiate = TRUE, conf.int = TRUE) %>% mutate(
  HR = sprintf("%.2f [%.2f, %.2f]", estimate, conf.low, conf.high)
)
# add risk difference and format model results
rd <- list()
for (x in ENTRY){ rd[[x]] <- cox_rd(cox_entry, survival_data, group_var = x, taus) }
rd <- bind_rows(rd)
cox.entry.tidy <- cox.entry.tidy %>% left_join(rd, by = "term") %>% rename(fHR = HR, HR = estimate) %>% 
  dplyr::select(-risk) %>% pivot_wider(names_from = tau, names_prefix = "RD_at_act", values_from = risk_difference) # pivot taus as columns
cox.entry.res <- cox.entry.tidy %>% dplyr::select(term, HR)
# bootstrap call: may take a few minutes
boots.entry <- cox_boot_ci(cox_entry_formula, survival_data, label = "Cox_entry")
boots.entry <- boots.entry %>% dplyr::select(-risk) %>% pivot_wider(names_from = tau, names_prefix = "RD_at_act", values_from = risk_difference) # pivot taus as columns
cox.entry.res <- add_boot_hr(cox.entry.res, cox.entry.tidy, boots.entry)
cox.entry.res

write.xlsx(cox.entry.res %>% as.data.frame(), file.path(models_dir, "models.xlsx"), sheetName = "Cox_entry", append = T, row.names = F, showNA = F)

#### cox extended adjustment ----
cox_extended_formula <- as.formula(paste0(
  "Surv(surv_time, surv_event) ~ ",
  COV_EXTENDED,
  " + cluster(play_title)"
))
cox_extended <- coxph(
  cox_extended_formula,
  data = survival_data,
  ties = "breslow"
)
survival::cox.zph(cox_extended)
cox.extended.tidy <- broom::tidy(cox_extended, exponentiate = TRUE, conf.int = TRUE) %>% mutate(
  HR = sprintf("%.2f [%.2f, %.2f]", estimate, conf.low, conf.high)
)
# add risk difference and format model results
rd <- list()
for (x in EXTENDED){ rd[[x]] <- cox_rd(cox_extended, survival_data, group_var = x, taus) }
rd <- bind_rows(rd)
cox.extended.tidy <- cox.extended.tidy %>% left_join(rd, by = "term") %>% rename(fHR = HR, HR = estimate) %>% 
  dplyr::select(-risk) %>% pivot_wider(names_from = tau, names_prefix = "RD_at_act", values_from = risk_difference) # pivot taus as columns
cox.extended.res <- cox.extended.tidy %>% dplyr::select(term, HR)
# bootstrap call: may take a few minutes
boots.extended <- cox_boot_ci(cox_extended_formula, survival_data, label = "Cox_extended")
boots.extended <- boots.extended %>% dplyr::select(-risk) %>% pivot_wider(names_from = tau, names_prefix = "RD_at_act", values_from = risk_difference) # pivot taus as columns
cox.extended.res <- add_boot_hr(cox.extended.res, cox.extended.tidy, boots.extended)
cox.extended.res

write.xlsx(cox.extended.res %>% as.data.frame(), file.path(models_dir, "models.xlsx"), sheetName = "Cox_extended", append = T, row.names = F, showNA = F)

#### HYPOTETHICAL: cox extended with time-varying supernatural exposure SIMULATED before ----
cox_extended_tv_formula <- as.formula(paste0(
  "Surv(surv_time, surv_event) ~ ",
  str_replace(COV_EXTENDED, "supernatural_exposure", " tt(tv_supernatural_exposure)"), # take supernatural exposure as time varying: tt(tv_supernatural_exposure)
  " + cluster(play_title)"
))
cox_extended_tv <- coxph(
  cox_extended_tv_formula,
  tt = tt_fun,
  data = survival_data,
  ties = "breslow"
)
summary(cox_extended_tv)
cox.extended.tidy.tv <- broom::tidy(cox_extended_tv, exponentiate = TRUE, conf.int = TRUE) %>% mutate(
  HR = sprintf("%.2f [%.2f, %.2f]", estimate, conf.low, conf.high)
)
write.xlsx(cox.extended.tidy.tv %>% as.data.frame(), file.path(models_dir, "models.xlsx"), sheetName = "Cox_extended_TV-supernatural", append = T, row.names = F, showNA = F)

## RMST  ----
tau_val <- 5  # time horizon to 5 acts
# function to run 2000 bootstrap samples RMST to build robust CIs
rmst_by_group_boot <- function(data, time_var, status_var, group_var, tau, adjust_vars = NULL, n_boot = 2000) {
  grp <- data[[group_var]] # take current variable (e.g. social status, genre, ...) 
  levs <- levels(as.factor(grp)) # take its levels
  if (length(levs) < 2) return(NULL) 
  ref <- levs[1]; others <- levs[-1] # separate reference level from others
  play_ids <- unique(data$play_title) # take unique play titles (cluster var)
  
  fit_one <- function(d, lev) { # original sample estimated RMST
    arm <- as.numeric(d[[group_var]] == lev) # each level of group_var at a time
    cov_vars <- setdiff(adjust_vars, group_var) # select adjustment vars (everything else, minus group_var itself)
    #cov_vars <- setdiff(cov_vars, "genre") # drop genre: linearly dependent with year_cut, kept out of the design matrix
    covs_df <- if (length(cov_vars) > 0) # build model matrix 
      model.matrix(~ ., data = d[, cov_vars, drop = FALSE])[, -1, drop = FALSE] else NULL
    tryCatch( # try Catch needed not to stop due to errors (bootstrap sampling may fail a model by chance ...)
      rmst2(time = d[[time_var]], status = d[[status_var]], arm = arm, tau = tau, covariates = covs_df), 
      # initialise RMST model with time variable (surv_time), status (surv_event), group_var level (arm), time horizon (tau), and adjustments (covs)
      error = function(e) NULL
    )
  }
  
  # pull every term's estimate out of an rmst2 fit
  extract_all <- function(f, adjusted) {
    if (adjusted) {
      dplyr::bind_rows( # stack the three RMST result matrices into one long term/stat/value table
        as.data.frame(f$RMST.difference.adjusted) %>% tibble::rownames_to_column("term") %>%
          dplyr::transmute(term, stat = "diff", value = coef), # difference in survival time, one row per term (arm + each adjust_var dummy)
        as.data.frame(f$RMST.ratio.adjusted) %>% tibble::rownames_to_column("term") %>%
          dplyr::transmute(term, stat = "rat", value = `exp(coef)`), # ratio of survival time, one row per term
        as.data.frame(f$RMTL.ratio.adjusted) %>% tibble::rownames_to_column("term") %>%
          dplyr::transmute(term, stat = "rmtl", value = `exp(coef)`) # ratio of time-lost, one row per term
      )
    } else {
      # if no-adjustment --> RMST gives a different results structure, and only "arm" exists (no adjust_vars to report)
      u <- f$unadjusted.result
      tibble::tibble(
        term  = "arm",
        stat  = c("diff", "rat", "rmtl"),
        value = c(u["RMST (arm=1)-(arm=0)", "Est."],
                  u["RMST (arm=1)/(arm=0)", "Est."],
                  u["RMTL (arm=1)/(arm=0)", "Est."])
      )
    }
  }
  
  purrr::map_dfr(others, function(lvl) { # do this - map - across all levels of group_var compared to ref
    sub_data <- data[grp %in% c(ref, lvl), ] # subset with this turn lvl vs ref
    adjusted <- !is.null(adjust_vars) && length(setdiff(setdiff(adjust_vars, group_var), "genre")) > 0 # flag: are there adjustment vars beyond group_var and genre
    
    point <- fit_one(sub_data, lvl) # point estimate
    if (is.null(point)) return(NULL) # if model fails, returns NULL
    point_est <- extract_all(point, adjusted) # point estimate for THIS group_var's arm term only (adjustment terms discarded downstream, see wrapper)
    
    set.seed(2604) # make random sampling reproducible (to expand, see https://www.uvm.edu/~cbcafier/cs1210/book/12_randomness/seed.html)
    boots <- purrr::map_dfr(seq_len(n_boot), function(i) { # may take a while depending on your PC characteristics, do RMST models on 2000 samples
      sampled_plays <- sample(play_ids, length(play_ids), replace = TRUE) # (pseudo-) randomly sample plays
      d_boot <- purrr::map_dfr(sampled_plays, ~ sub_data[sub_data$play_title == .x, ]) # b-th sample dataset
      f <- fit_one(d_boot, lvl) # just fit the model on this dataset now
      if (is.null(f)) return(NULL) #return NULL if failing model
      extract_all(f, adjusted) %>% dplyr::mutate(iteration = i) # same extraction as point estimate, tagged with iteration number
    }, .progress = paste("Bootstrap: ", group_var, "-", lvl)) # progress bar for each lvl
    # boots is the bootstrap distribution of diff, rat and rmtl for EVERY term: i.e. 2,000 possible versions of your model - they let you guess to what extent your findings may vary by random chance
    
    quant <- function(x, p) stats::quantile(x, p, na.rm = TRUE, names = FALSE) # take boot distribution quantiles
    boots_ci <- boots %>%
      dplyr::group_by(term, stat) %>% # compute CI separately per term (arm, each adjust_var) and per stat (diff/rat/rmtl)
      dplyr::summarise(lo = quant(value, .025), hi = quant(value, .975), .groups = "drop")
    
    point_est %>%
      dplyr::left_join(boots_ci, by = c("term", "stat")) %>% # attach bootstrap CI to each point estimate
      dplyr::mutate(est_ci = sprintf("%.2f [%.2f, %.2f]", value, lo, hi)) %>% # format as estimate [95%CI]
      dplyr::select(term, stat, est_ci) %>%
      tidyr::pivot_wider(names_from = stat, values_from = est_ci,
                         names_glue = "RMST{toupper(substr(stat,1,1))}{substr(stat,2,nchar(stat))}") %>% # reshape diff/rat/rmtl into RMSTDiff/RMSTRat/RMSTRmtl columns
      dplyr::mutate(Variable = group_var, group_level = lvl, # exported results
                    n_boot_ok = dplyr::n_distinct(boots$iteration)) %>% # how many models did not fail (see above)
      dplyr::relocate(Variable, group_level, term) # just move columns in desired order
  })
}

# function to calculate bootstrap RMST CI for each variable adjusted by the others
rmst_all_vars_boot <- function(data, time_var, status_var, vars_to_test, tau, n_boot = 2000) {
  purrr::map_dfr(vars_to_test, function(v) { # loop: each variable takes a turn as group_var
    adj <- setdiff(vars_to_test, v) # every other variable in the list becomes its adjustment set
    out <- rmst_by_group_boot(
      data = data, time_var = time_var, status_var = status_var,
      group_var = v, adjust_vars = adj, tau = tau, n_boot = n_boot
    )
    if (is.null(out)) return(NULL)
    out %>% dplyr::filter(term == "arm") # keep only the group_var's own coefficient (its adjustment-term rows come from a model where IT was the arm, not from THEIR own turn as group_var, so they're discarded here to avoid duplicated/mismatched rows)
  })
}

# now use previous function to estimate RMST with 2,000-samples bootstrap 95%CI 
rmst_results_unadj <- rmst_all_vars_boot(survival_data, time_var = "surv_time", status_var = "surv_event", tau = tau_val,
                                         vars_to_test = "social_status", # crude model
                                         n_boot = 2000)
rmst_results_entryadj <- rmst_all_vars_boot(survival_data, time_var = "surv_time", status_var = "surv_event", tau = tau_val,
                                            vars_to_test = c("social_status", "genre", "gender", "year_cut"), #entry adjustment 
  # NOTE: year_cut is linearly dependent with genre and cannot be handled properly by RMST function as it is by coxph since its Newton-Raphson estimation is fairly tolerant of near-singular columns
                                            n_boot = 2000)
rmst_results_extendedadj <- rmst_all_vars_boot(survival_data, time_var = "surv_time", status_var = "surv_event", tau = tau_val,
                                               vars_to_test = c("social_status", "genre", "gender", "year_cut", "character_status", "supernatural_exposure", "betrayal_exposure_cat"), # extended adjustment
                                               # NOTE: year_cut is linearly dependent with genre and cannot be handled properly by RMST function as it is by coxph since its Newton-Raphson estimation is fairly tolerant of near-singular columns
                                               n_boot = 2000)
# specify which model findings
rmst_results_unadj$model <- "Unadjusted RMST"
rmst_results_entryadj$model <- "Entry-adjusted RMST"
rmst_results_extendedadj$model <- "Extended-adjusted RMST"
# merge together
rmst_results <- bind_rows(rmst_results_unadj, rmst_results_entryadj, rmst_results_extendedadj) %>% as.data.frame()
# write results on models.xlsx
write.xlsx(rmst_results, file = file.path(models_dir, "models.xlsx"), sheetName = "RMST", append = T, row.names = F, showNA = F)

## Pooled logistic model: standardized cumulative mortality ----
pooled_base <- survival_data %>%
  dplyr::select(
    character_id, play_title, surv_time, surv_event,
    social_status, genre, gender, character_status,
    supernatural_exposure, betrayal_exposure_cat, year_cut
  ) %>%
  dplyr::filter(
    dplyr::if_all(dplyr::everything(), ~ !is.na(.))
  ) %>%
  dplyr::mutate(
    # Assign observations to one of the five discrete character-act intervals
    followup_act = pmin(5L, pmax(1L, ceiling(as.numeric(surv_time)))),
    social_status = relevel_if_present(social_status, "Commoner")
  )
status_levels <- levels(pooled_base$social_status)
# Convert character-level data into character-act data
expand_character_acts <- function(data) {
  data %>%
    tidyr::uncount(
      weights = followup_act,
      .remove = FALSE,
      .id = "act_number"
    ) %>%
    dplyr::mutate(
      interval_event = as.integer(
        surv_event == 1 & act_number == followup_act
      ),
      act_factor = factor(act_number, levels = 1:5)
    )
}
pooled_data <- expand_character_acts(pooled_base)
pooled_early <- pooled_data %>% dplyr::filter(act_number <= 3)
pooled_late  <- pooled_data %>% dplyr::filter(act_number >= 4)

# Note: droplevels() so act_factor has only the levels present in each piece
pooled_early <- pooled_early %>% dplyr::mutate(act_factor = droplevels(act_factor))
pooled_late  <- pooled_late  %>% dplyr::mutate(act_factor = droplevels(act_factor))


#### Helpers: fit the two pieces and combine into one curve ----
# set: set of adjusting covariates, entry or extended model
fit_pooled <- function(d, set) {
  glm(
    as.formula(paste0("interval_event ~ ", paste0(set, collapse = " + "))),
    data = d, family = binomial("logit")
  )
}

tidy_cl <- function(fit, d) {
  vc <- sandwich::vcovCL(fit, cluster = d$play_title)
  lmtest::coeftest(fit, vcov = vc) %>%
    broom::tidy(conf.int = TRUE) %>%
    dplyr::mutate(OR = sprintf("%.2f [%.2f, %.2f]", exp(estimate), exp(conf.low), exp(conf.high)), p.value = round(p.value, 3)) %>%
    dplyr::select(term, OR, p.value)
}
fit_split <- function(long_data, set) {
  early <- long_data %>% dplyr::filter(act_number <= 3) %>%
    dplyr::mutate(act_factor = droplevels(act_factor))
  late  <- long_data %>% dplyr::filter(act_number >= 4) %>%
    dplyr::mutate(act_factor = droplevels(act_factor))
  list(early = fit_pooled(early, set), late = fit_pooled(late, set))
}

standardized_curve_split <- function(fit_e, fit_l, data) {
  s_e <- standardized_mortality(fit_e, data, acts = 1:3)
  s_l <- standardized_mortality(fit_l, data, acts = 4:5) %>%
    dplyr::left_join(
      s_e %>% dplyr::filter(act_number == 3) %>%
        dplyr::select(standardization_id, level, s3 = surv),
      by = c("standardization_id", "level")
    ) %>%
    dplyr::mutate(surv = surv * s3) %>%   # S(t) = S_early(3) * S_late(t | alive at 3)
    dplyr::select(-s3)
  dplyr::bind_rows(s_e, s_l) %>%
    dplyr::group_by(social_status = level, act_number) %>%
    dplyr::summarise(cumulative_mortality = 1 - mean(surv), .groups = "drop")
}


bootstrap_once <- function(iteration, set) {
  sampled_plays <- sample(play_ids, size = length(play_ids), replace = TRUE) 
  
  idx_list <- idx_by_play[as.character(sampled_plays)]
  bootstrap_number <- rep(seq_along(sampled_plays), lengths(idx_list))
  bootstrap_data <- pooled_base[unlist(idx_list, use.names = FALSE), ]
  bootstrap_data$bootstrap_play <- bootstrap_number
  bootstrap_data$character_id <- paste0(bootstrap_number, "_", bootstrap_data$character_id)
  
  bootstrap_long <- expand_character_acts(bootstrap_data)
  
  bootstrap_curve <- tryCatch(
    suppressWarnings({
      bf <- fit_split(bootstrap_long, set)
      standardized_curve_split(bf$early, bf$late, bootstrap_data)
    }),
    error = function(e) NULL
  )
  if (is.null(bootstrap_curve)) return(NULL)
  ref <- bootstrap_curve[bootstrap_curve$social_status == "Commoner", ]
  bootstrap_curve$commoner_mortality <- ref$cumulative_mortality[match(bootstrap_curve$act_number, ref$act_number)]
  bootstrap_curve$risk_difference <- bootstrap_curve$cumulative_mortality - bootstrap_curve$commoner_mortality
  bootstrap_curve$iteration <- iteration
  bootstrap_curve
}
quantile_safe <- function(x, probability) {
  if (all(is.na(x))) {
    return(NA_real_)
  }
  stats::quantile(
    x,
    probs = probability,
    na.rm = TRUE,
    names = FALSE
  )
}

#### -----

### fit crude ----
fit_early <- fit_pooled(pooled_early, CRUDE)
fit_late  <- fit_pooled(pooled_late, CRUDE)

res_early <- tidy_cl(fit_early, pooled_early) %>% mutate(acts = "1-3")
res_late  <- tidy_cl(fit_late,  pooled_late) %>% mutate(acts = "4-5")

pooledlogit.crude.res <- bind_rows(res_early, res_late)
write.xlsx(pooledlogit.crude.res %>% as.data.frame() %>% mutate(p.value = round(p.value, 3)), file.path(models_dir, "models.xlsx"), sheetName = "Pooled_logistic_crude", append = T, row.names = F, showNA = F)

#### Standardized cumulative mortality curves ----
# to obtain absolute risk measures from pooled logistic regression
standardized_mortality <- function(model, character_data, acts, var = "social_status", levels_ = status_levels) {
  purrr::map_dfr(levels_, function(v) {
    nd <- character_data %>%
      dplyr::mutate(standardization_id = dplyr::row_number()) %>%
      tidyr::uncount(weights = length(acts), .remove = FALSE) %>%
      dplyr::group_by(standardization_id) %>%
      dplyr::mutate(act_number = acts) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        act_factor = factor(act_number, levels = acts),
        !!var := if (is.factor(character_data[[var]])) factor(v, levels = levels_) else v
      )
    nd$h <- predict(model, newdata = as.data.frame(nd), type = "response")
    nd %>%
      dplyr::arrange(standardization_id, act_number) %>%
      dplyr::group_by(standardization_id) %>%
      dplyr::mutate(surv = cumprod(1 - h)) %>%   # conditional survival within the piece
      dplyr::ungroup() %>%
      dplyr::select(standardization_id, act_number, surv) %>%
      dplyr::mutate(level = v)
  })
}

s_early <- standardized_mortality(fit_early, pooled_base, acts = 1:3)
s_late  <- standardized_mortality(fit_late,  pooled_base, acts = 4:5)

# combine per individual: S(t) = S_early(3) * S_late(t), then average
res <- dplyr::bind_rows(
  s_early,
  s_late %>% dplyr::left_join(
    s_early %>% dplyr::filter(act_number == 3) %>%
      dplyr::select(standardization_id, level, s3 = surv),
    by = c("standardization_id", "level")) %>%
    dplyr::mutate(surv = surv * s3) %>% dplyr::select(-s3)
) %>%
  dplyr::group_by(level, act_number) %>%
  dplyr::summarise(cum_mortality = 1 - mean(surv), .groups = "drop")



# Point estimate
fits <- fit_split(pooled_data, CRUDE)
standardized_curve <- standardized_curve_split(fits$early, fits$late, pooled_base) %>%
  dplyr::group_by(act_number) %>%
  dplyr::mutate(
    commoner_mortality = cumulative_mortality[social_status == "Commoner"][1],
    risk_difference    = cumulative_mortality - commoner_mortality
  ) %>%
  dplyr::ungroup()

# Cluster bootstrap by play 
set.seed(20260927)
n_bootstrap <- 2000
play_ids <- unique(pooled_base$play_title)

idx_by_play <- base::split(seq_len(nrow(pooled_base)), pooled_base$play_title)[as.character(play_ids)]
bootstrap_results <- purrr::map_dfr(seq_len(n_bootstrap), bootstrap_once, set = CRUDE, .progress = "b: ")

# Bootstrap confidence intervals 
bootstrap_ci <- bootstrap_results %>%
  dplyr::group_by(social_status, act_number) %>%
  dplyr::summarise(
    mortality_lower =
      quantile_safe(cumulative_mortality, 0.025),
    mortality_upper =
      quantile_safe(cumulative_mortality, 0.975),
    difference_lower =
      quantile_safe(risk_difference, 0.025),
    
    difference_upper =
      quantile_safe(risk_difference, 0.975),
    .groups = "drop"
  )
absolute_risk_results <- standardized_curve %>%
  dplyr::left_join(
    bootstrap_ci,
    by = c("social_status", "act_number")
  )
# Supplementary table: estimates at three and five character-acts
supplementary_table_absolute_risk <- absolute_risk_results %>%
  dplyr::filter(act_number %in% c(3, 5)) %>%
  dplyr::mutate(
    `Cumulative mortality x100, % [95% CI]` = sprintf(
      "%.1f [%.1f, %.1f]",
      100 * cumulative_mortality,
      100 * mortality_lower,
      100 * mortality_upper
    ),
    `Absolute risk difference x100, percentage points [95% CI]` =
      dplyr::if_else(
        social_status == "Commoner",
        "Reference",
        sprintf(
          "%+.1f [%.1f, %.1f]",
          100 * risk_difference,
          100 * difference_lower,
          100 * difference_upper
        )
      )
  ) %>%
  dplyr::select(
    `Social status` = social_status,
    `Character-act` = act_number,
    `Cumulative mortality x100, % [95% CI]`,
    `Absolute risk difference x100, percentage points [95% CI]`
  ) %>%
  dplyr::arrange(`Character-act`, `Social status`)
supplementary_table_absolute_risk
# save table 
xlsx::write.xlsx(
  as.data.frame(supplementary_table_absolute_risk), file = file.path( models_dir, "models.xlsx"),
  sheetName = "Absolute risk table_pl_crude", row.names = FALSE, append = TRUE, showNA = FALSE
)

### fit entry adjustment model ----
# Point estimate
fits <- fit_split(pooled_data, COV_ENTRY)
standardized_curve <- standardized_curve_split(fits$early, fits$late, pooled_base) %>%
  dplyr::group_by(act_number) %>%
  dplyr::mutate(
    commoner_mortality = cumulative_mortality[social_status == "Commoner"][1],
    risk_difference    = cumulative_mortality - commoner_mortality
  ) %>%
  dplyr::ungroup()

# Cluster bootstrap by play 
set.seed(20260927)
n_bootstrap <- 2000
play_ids <- unique(pooled_base$play_title)

idx_by_play <- base::split(seq_len(nrow(pooled_base)), pooled_base$play_title)[as.character(play_ids)]
bootstrap_results <- purrr::map_dfr(seq_len(n_bootstrap), bootstrap_once, set = COV_ENTRY, .progress = "b: ")

# Bootstrap confidence intervals 
bootstrap_ci <- bootstrap_results %>%
  dplyr::group_by(social_status, act_number) %>%
  dplyr::summarise(
    mortality_lower =
      quantile_safe(cumulative_mortality, 0.025),
    mortality_upper =
      quantile_safe(cumulative_mortality, 0.975),
    difference_lower =
      quantile_safe(risk_difference, 0.025),
    
    difference_upper =
      quantile_safe(risk_difference, 0.975),
    .groups = "drop"
  )
absolute_risk_results <- standardized_curve %>%
  dplyr::left_join(
    bootstrap_ci,
    by = c("social_status", "act_number")
  )
# Supplementary table: estimates at three and five character-acts
supplementary_table_absolute_risk <- absolute_risk_results %>%
  dplyr::filter(act_number %in% c(3, 5)) %>%
  dplyr::mutate(
    `Cumulative mortality x100, % [95% CI]` = sprintf(
      "%.1f [%.1f, %.1f]",
      100 * cumulative_mortality,
      100 * mortality_lower,
      100 * mortality_upper
    ),
    `Absolute risk difference x100, percentage points [95% CI]` =
      dplyr::if_else(
        social_status == "Commoner",
        "Reference",
        sprintf(
          "%+.1f [%.1f, %.1f]",
          100 * risk_difference,
          100 * difference_lower,
          100 * difference_upper
        )
      )
  ) %>%
  dplyr::select(
    `Social status` = social_status,
    `Character-act` = act_number,
    `Cumulative mortality x100, % [95% CI]`,
    `Absolute risk difference x100, percentage points [95% CI]`
  ) %>%
  dplyr::arrange(`Character-act`, `Social status`)
supplementary_table_absolute_risk
# save table 
xlsx::write.xlsx(
  as.data.frame(supplementary_table_absolute_risk), file = file.path( models_dir, "models.xlsx"),
  sheetName = "Absolute risk table_pl_entry", row.names = FALSE, append = TRUE, showNA = FALSE
)

### fit extended adjustment ----
fit_early <- fit_pooled(pooled_early, COV_EXTENDED)
fit_late  <- fit_pooled(pooled_late, COV_EXTENDED)

res_early <- tidy_cl(fit_early, pooled_early) %>% mutate(acts = "1-3")
res_late  <- tidy_cl(fit_late,  pooled_late) %>% mutate(acts = "4-5")

pooledlogit.ext.res <- bind_rows(res_early, res_late)
write.xlsx(pooledlogit.ext.res %>% as.data.frame() %>% mutate(p.value = round(p.value, 3)), file.path(models_dir, "models.xlsx"), sheetName = "Pooled_logistic_ext_adj", append = T, row.names = F, showNA = F)

# Point estimate
fits <- fit_split(pooled_data, COV_EXTENDED)
standardized_curve <- standardized_curve_split(fits$early, fits$late, pooled_base) %>%
  dplyr::group_by(act_number) %>%
  dplyr::mutate(
    commoner_mortality = cumulative_mortality[social_status == "Commoner"][1],
    risk_difference    = cumulative_mortality - commoner_mortality
  ) %>%
  dplyr::ungroup()

# Cluster bootstrap by play 
set.seed(20260927)
n_bootstrap <- 2000
play_ids <- unique(pooled_base$play_title)

bootstrap_results <- purrr::map_dfr(seq_len(n_bootstrap), bootstrap_once, set = COV_EXTENDED, .progress = "b: ")

# Bootstrap confidence intervals
bootstrap_ci <- bootstrap_results %>%
  dplyr::group_by(social_status, act_number) %>%
  dplyr::summarise(
    mortality_lower =
      quantile_safe(cumulative_mortality, 0.025),
    mortality_upper =
      quantile_safe(cumulative_mortality, 0.975),
    difference_lower =
      quantile_safe(risk_difference, 0.025),
    
    difference_upper =
      quantile_safe(risk_difference, 0.975),
    .groups = "drop"
  )
absolute_risk_results <- standardized_curve %>%
  dplyr::left_join(
    bootstrap_ci,
    by = c("social_status", "act_number")
  )
# Supplementary table: estimates at three and five character-acts
supplementary_table_absolute_risk <- absolute_risk_results %>%
  dplyr::filter(act_number %in% c(3, 5)) %>%
  dplyr::mutate(
    `Cumulative mortality x100, % [95% CI]` = sprintf(
      "%.1f [%.1f, %.1f]",
      100 * cumulative_mortality,
      100 * mortality_lower,
      100 * mortality_upper
    ),
    `Absolute risk difference x100, percentage points [95% CI]` =
      dplyr::if_else(
        social_status == "Commoner",
        "Reference",
        sprintf(
          "%+.1f [%.1f, %.1f]",
          100 * risk_difference,
          100 * difference_lower,
          100 * difference_upper
        )
      )
  ) %>%
  dplyr::select(
    `Social status` = social_status,
    `Character-act` = act_number,
    `Cumulative mortality x100, % [95% CI]`,
    `Absolute risk difference x100, percentage points [95% CI]`
  ) %>%
  dplyr::arrange(`Character-act`, `Social status`)
supplementary_table_absolute_risk
# save table 
xlsx::write.xlsx(
  as.data.frame(supplementary_table_absolute_risk), file = file.path( models_dir, "models.xlsx"),
  sheetName = "Absolute risk table_pl_ext", row.names = FALSE, append = TRUE, showNA = FALSE
)

# Supplementary figure ----
plot_data <- absolute_risk_results %>% # takes the last model: extended adjustment
  dplyr::select(
    social_status, act_number,
    cumulative_mortality,
    mortality_lower, mortality_upper
  )
# Add time zero
plot_data <- dplyr::bind_rows(
  tibble::tibble(
    social_status = status_levels,
    act_number = 0L,
    cumulative_mortality = 0,
    mortality_lower = 0,
    mortality_upper = 0
  ),
  plot_data
) %>%
  dplyr::mutate(
    social_status = factor(
      social_status,
      levels = status_levels
    )
  )
supplementary_figure_absolute_risk <- ggplot2::ggplot(
  plot_data,
  ggplot2::aes(
    x = act_number,
    y = cumulative_mortality,
    colour = social_status,
    group = social_status
  )
) +
  ggplot2::geom_ribbon(
    ggplot2::aes(
      ymin = mortality_lower,
      ymax = mortality_upper,
      fill = social_status
    ),
    alpha = 0.15,
    colour = NA
  ) +
  ggplot2::geom_step(
    linewidth = 1
  ) +
  ggplot2::scale_x_continuous(
    breaks = 0:5,
    limits = c(0, 5)
  ) +
  ggplot2::scale_y_continuous(
    breaks = seq(0, .5, by = 0.1),
    limits = c(0, 1)
  ) +
  ggplot2::labs(
    x = "Character-act",
    y = "Standardized cumulative mortality",
    colour = "Social status",
    fill = "Social status"
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    legend.position = "bottom",
    panel.grid.minor = ggplot2::element_blank()
  ) + facet_wrap(~ social_status)
supplementary_figure_absolute_risk
ggplot2::ggsave(
  filename = file.path(
    figures_dir,
    "Supplementary_standardized_cumulative_mortality.png"
  ),
  plot = supplementary_figure_absolute_risk,
  width = 8,
  height = 6,
  dpi = 300
)
# Figures ----
## bar plot for mortality rates by play, genre, by profession ----
fig_play_rates <- by_play %>%
  mutate(play_title = forcats::fct_reorder(play_title, incidence_rate_per_100_character_acts)) %>%
  ggplot(aes(x = play_title, y = incidence_rate_per_100_character_acts, fill = genre)) +
  geom_col() +
  scale_fill_manual(values = col_genre, name = "Genre")+
  coord_flip() +
  labs(
    title = "Mortality rate by play",
    x = NULL,
    y = "Deaths per 100 character-acts"
  ) + geom_errorbar(
    aes(ymin = incidence_rate_lower_95, ymax = incidence_rate_upper_95),
    width = 0.2
  ) +
  theme_minimal(base_size = 12) + 
  theme(plot.title = element_text(family = font2, size = 36, hjust = .5), 
        axis.text.y = element_text(family = font2, size = 30, face = "plain"),
        axis.title.x = element_text(family = font2, size = 24),
        axis.text.x = element_text(family = font2, size = 24, face = "plain"),
        legend.title    = element_text(family = font2, size = 24),
        legend.text     = element_text(family = font2, size = 20), legend.position = "none")
fig_play_rates

# Incidence rate by profession group
fig_profession_rates <- by_profession_group %>% filter(profession_group != "Other or unclassified") %>%
  mutate(profession_group = forcats::fct_reorder(profession_group, incidence_rate_per_100_character_acts)) %>%
  ggplot(aes(x = profession_group, y = incidence_rate_per_100_character_acts, fill = profession_group)) +
  geom_col() +
  scale_fill_manual(values = paste0("grey", rev(seq(10, 90, by = 10))), name = "Profession")+
  geom_errorbar(
    aes(ymin = incidence_rate_lower_95, ymax = incidence_rate_upper_95),
    width = 0.2
  ) +
  coord_flip() +
  labs(
    title = "Mortality rate by profession group",
    x = NULL,
    y = "Deaths per 100 character-acts"
  ) +
  theme_minimal(base_size = 12)+
  theme(plot.title = element_text(family = font2, size = 36, hjust = .5), 
        axis.text.y = element_text(family = font2, size = 30, face = "plain"),
        axis.title.x = element_text(family = font2, size = 24),
        axis.text.x = element_text(family = font2, size = 24, face = "plain"),
        legend.title    = element_text(family = font2, size = 24),
        legend.text     = element_text(family = font2, size = 20), legend.position = "none")
fig_profession_rates

# incidence rate by genre
fig_genre_rates <- by_genre %>% 
  mutate(genre = forcats::fct_reorder(genre, incidence_rate_per_100_character_acts)) %>%
  ggplot(aes(x = genre, y = incidence_rate_per_100_character_acts, fill = genre)) +
  geom_col() +
  scale_fill_manual(values = col_genre, name = "Genre")+
  geom_errorbar(
    aes(ymin = incidence_rate_lower_95, ymax = incidence_rate_upper_95),
    width = 0.2
  ) +
  coord_flip() +
  labs(
    title = "Mortality rate by genre",
    x = NULL,
    y = "Deaths per 100 character-acts"
  ) +
  theme_minimal(base_size = 12)+
  theme(plot.title = element_text(family = font2, size = 36, hjust = .5), 
        axis.text.y = element_text(family = font2, size = 30, face = "plain"),
        axis.title.x = element_text(family = font2, size = 24),
        axis.text.x = element_text(family = font2, size = 24, face = "plain"),
        legend.title    = element_text(family = font2, size = 24),
        legend.text     = element_text(family = font2, size = 20), legend.position = "none")
fig_genre_rates

bg_theme <- theme(
  plot.background  = element_rect(fill = bg_color, color = NA),
  panel.background = element_rect(fill = bg_color, color = NA),
  legend.background = element_rect(fill = bg_color, color = NA),
  legend.box.background = element_rect(fill = bg_color, color = NA)
)
FIG1.rates <- fig_play_rates / (fig_profession_rates + fig_genre_rates) + plot_layout(heights = c(2, 1)) +
  plot_annotation(theme = theme(plot.background = element_rect(fill = bg_color, color = NA))) & bg_theme
FIG1.rates
jpeg(file = paste0(figures_dir, "/RATES.jpg"), width = 10, height = 7, units = "in", res = 300, quality = 100)
print(FIG1.rates); dev.off()

## map ----
library(rnaturalearth); library(sf); library(dplyr); library(tidyr); library(ggplot2); library(ggrepel)

world <- ne_countries(scale = "medium", returnclass = "sf")
# north africa and middle east
north_africa_me <- c("Algeria", "Armenia", "Azerbaijan", "Chad", "Cyprus", "Egypt", "Georgia", "Iran", "Iraq", "Israel", "Jordan", "Lebanon", "Libya", 
                     "Mali", "Mauritania", "Morocco", "Niger", "Saudi Arabia", "Senegal", "Sudan", "Syria", "Tunisia", "Turkey", "Western Sahara")
europe_na <- world %>%
  filter(continent == "Europe" | name %in% north_africa_me) %>%
  filter(name != "Russia" | TRUE)
by_country <- analysis_data %>%
  mutate(country = ifelse(play_title %in% "The Winter's Tale", "Czechia", country),
         country = ifelse(play_title %in% "The Tempest", "Italy", country)
  ) %>% group_by(country) %>% summarise(
    n_characters = n(),
    play_titles = paste(unique(play_title), collapse = "\n"),
    n_plays = n_distinct(play_title),
    character_acts_total = sum(character_acts, na.rm = TRUE),
    deaths = sum(death_event, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  add_rate_columns() %>%
  arrange(desc(deaths), desc(incidence_rate_per_100_character_acts), country)
europe_na <- europe_na %>% left_join(by_country %>% 
                                       rename(mortality_rate = incidence_rate_per_100_character_acts), by = join_by(name_ciawf == country))
europe_na$name_ciawf[!is.na(europe_na$mortality_rate)]
base_map <- ggplot(europe_na) +
  geom_sf(aes(fill = mortality_rate), color = "grey80", alpha = .6, linewidth = 0.25) +
  coord_sf(xlim = c(-25, 45), ylim = c(27, 60), expand = FALSE) +
  scale_fill_gradientn(
    colours = c("#003559", "red4"),
    na.value = "white",
    name = "Deaths per\n100 character-acts"
  ) +
  theme_void(base_size = 12)
base_map

col_genre <- c("History" = "orange", "Romance" = "black", "Tragedy" = "red4")
# --- precise settings (approximate coordinates for the play's actual location) ---
play_settings <- tribble(
  ~play_title,                    ~lon,      ~lat,
  "Antony and Cleopatra",         29.9187,   31.2001,  # Alexandria, Egypt
  "Coriolanus",                   12.4964,   41.9028,  # Rome, Italy
  "Cymbeline",                    -5.0341,   51.7128,  # Milford Haven, Wales
  "Hamlet",                       12.6091,   56.0399,  # Elsinore, Denmark
  "Henry IV, Part 1",             -2.7115,   52.7069,  # Shrewsbury (decisive battle)
  "Henry IV, Part 2",             -0.1278,   51.5074,  # London
  "Henry V",                       2.1275,   50.4547,  # Agincourt, France
  "Henry VI, Part 1",              1.9039,   47.9029,  # Orl\u00e9ans, France (siege)
  "Henry VI, Part 2",              0.0000,   52.5000,  # St Albans (1st battle)
  "Henry VI, Part 3",             -1.3336,   53.8266,  # Towton, Yorkshire
  "Henry VIII",                   -0.1278,   51.5074,  # London
  "Julius Caesar",                12.4964,   41.9028,  # Rome, Italy
  "King John",                    -0.5632,   47.4784,  # Angers, France
  "King Lear",                    -1.1398,   52.6369,  # Leicester (legendary seat)
  "Macbeth",                      -4.2247,   57.4791,  # Inverness, Scotland
  "Othello",                      33.4299,   35.1264,  # Cyprus
  "Pericles, Prince of Tyre",     35.2038,   33.2704,  # Tyre, Lebanon
  "Richard II",                   -0.1278,   51.5074,  # London
  "Richard III",                  -1.4058,   52.5813,  # Bosworth Field
  "Romeo and Juliet",             10.9916,   45.4384,  # Verona, Italy
  "The Tempest",                  14.2452,   37.0685,   # unnamed mediterranean island is unknown, assign Sicily
  "The Two Noble Kinsmen",        23.7275,   37.9838,  # Athens, Greece
  "The Winter's Tale",            14.4378,   50.0755,  # "Bohemia" (Prague; Shakespeare's famous coastal error)
  "Timon of Athens",              23.7275,   37.9838,  # Athens, Greece
  "Titus Andronicus",             12.4964,   41.9028,  # Rome, Italy
  "Troilus and Cressida",         26.2389,   39.9573,  # Troy, Turkey
)

# --- 1. unpack one row per play (as before) ---
play_geo <- europe_na %>%
  st_drop_geometry() %>%
  dplyr::select(name_ciawf, play_titles) %>%
  filter(!is.na(play_titles)) %>%
  tidyr::separate_rows(play_titles, sep = "\n") %>%
  rename(country = name_ciawf, play_title = play_titles) %>%
  left_join(europe_na %>% dplyr::select(name_ciawf, geometry), by = c("country" = "name_ciawf")) %>%
  st_as_sf()

# --- 2. jittered fallback point (for any play NOT in play_settings) ---
set.seed(42)

jitter_within <- function(geom, n) {
  pts <- sf::st_sample(sf::st_union(geom), size = n, type = "random", exact = TRUE)
  if (length(pts) < n) {
    extra <- sf::st_sample(sf::st_union(geom), size = n - length(pts), type = "random")
    pts <- c(pts, extra)
  }
  sf::st_coordinates(pts)
}

fallback_labels <- play_geo %>%
  group_by(country) %>%
  group_modify(~ {
    coords <- jitter_within(.x$geometry, nrow(.x))
    .x %>% st_drop_geometry() %>% mutate(lon_fallback = coords[, 1], lat_fallback = coords[, 2])
  }) %>%
  ungroup()

# --- 3. precise location wins where known, jitter only fills the gaps ---
play_labels <- fallback_labels %>%
  left_join(play_settings, by = "play_title") %>%
  mutate(
    lon = coalesce(lon, lon_fallback),
    lat = coalesce(lat, lat_fallback)
  ) %>%
  dplyr::select(country, play_title, lon, lat) %>%
  distinct(play_title, .keep_all = TRUE) %>%
  left_join(raw %>% dplyr::select(play_title, genre), by = "play_title") %>% distinct()

# --- 4. plot ---
fig_map <- base_map +
  geom_point(
    data = play_labels,
    aes(x = lon, y = lat, shape = genre, color = genre),
    inherit.aes = FALSE,
    size = 1.2,
    stroke = 0.8
  ) +
  scale_shape_manual(values = c("Tragedy" = 17, "History" = 15, "Romance" = 19), name = "Genre") +
  scale_color_manual(values = col_genre, name = "Genre") +
  geom_text_repel(
    data = play_labels,
    aes(x = lon, y = lat, label = play_title),
    inherit.aes = FALSE,
    family = "medieval",
    size = 7.5,
    color = "#2b2620",
    fontface = "bold",
    segment.color = "white",
    segment.size = 0.3,
    segment.alpha = 0.6,
    min.segment.length = 0,
    box.padding = 0.4,
    point.padding = 0.2,
    max.overlaps = Inf,
    seed = 42
  )

lon_range <- 45 - (-25)  
lat_range <- 60 - 27
mean_lat  <- (60 + 27) / 2

aspect_ratio <- (lon_range * cos(mean_lat * pi / 180)) / lat_range

fig_map <- fig_map +
  theme(
    plot.margin        = margin(t = 0, r = 50, b = 0, l = 10),
    legend.title        = element_text(family = "medieval", size = 18, lineheight = 0.55),
    legend.text         = element_text(family = "medieval", size = 18),
    plot.background     = element_rect(fill = "white"),
    panel.background    = element_rect(fill = "white"),
    legend.background   = element_rect(fill = "white")
  )
fig_map
jpeg(filename = paste0(figures_dir, "/MAP.jpg"), width = 7, height = 7 / aspect_ratio, units = "in", res = 300, quality = 100)
print(fig_map); dev.off()
# tables ----
table_1 <- analysis_data %>% summarise(
  n_characters = n(),
  n_plays = n_distinct(play_title),
  character_acts_total = sum(character_acts, na.rm = TRUE),
  deaths = sum(death_event, na.rm = TRUE),
  .groups = "drop"
) %>% add_rate_columns() %>%
  arrange(desc(deaths), desc(incidence_rate_per_100_character_acts))
table_1 <- data.frame(group_var = "overall", group_value = "overall", table_1)
group_vars <- c("genre", "social_status", "gender", "profession_group", "year_cut",
                "betrayal_exposure_cat", "supernatural_exposure",
                "character_status", "battle_related", "play_title")

summarise_by_var <- function(var, data) {
  data %>%
    group_by(group_value = .data[[var]]) %>%
    summarise(
      n_characters          = n(),
      n_plays               = n_distinct(play_title),
      character_acts_total  = sum(character_acts, na.rm = TRUE),
      deaths                = sum(death_event, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    add_rate_columns() %>%
    mutate(group_var = var, .before = group_value) %>%
    arrange(desc(incidence_rate_per_100_character_acts))
}

by_all_groups <- map_dfr(group_vars, summarise_by_var, data = analysis_data)

TABLE_1 <- bind_rows(table_1, by_all_groups) %>% mutate(
  MR = sprintf("%.2f [%.2f–%.2f]", incidence_rate_per_100_character_acts, incidence_rate_lower_95, incidence_rate_upper_95),
  deaths = sprintf("%.0f (%.1f%%)", deaths, cfr_percent)) %>%
  dplyr::select(-incidence_rate_per_100_character_acts, -incidence_rate_lower_95, -incidence_rate_upper_95, -cfr_percent)

write.xlsx(TABLE_1, file = paste0(tables_dir, "/TABLE_1.xlsx"), row.names = F, showNA = F)
