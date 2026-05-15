# Climate index preprocessing: download and merge NAO, AO, AMO, ENSO, PDO
setwd("data")

library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)

# --- Read NOAA-PSL monthly wide tables (year x 12 months) ---
read_psl_12 <- function(path, value_name) {
  lines <- readr::read_lines(path)
  lines <- lines[stringr::str_detect(lines, "^\\s*\\d{4}\\b")][-1]

  x <- readr::read_table(
    paste(lines, collapse = "\n"),
    col_names = c("year", month.abb),
    col_types = readr::cols(
      year = readr::col_integer(),
      Jan = readr::col_double(), Feb = readr::col_double(), Mar = readr::col_double(),
      Apr = readr::col_double(), May = readr::col_double(), Jun = readr::col_double(),
      Jul = readr::col_double(), Aug = readr::col_double(), Sep = readr::col_double(),
      Oct = readr::col_double(), Nov = readr::col_double(), Dec = readr::col_double()
    ),
    na = c("-99.99", "-99.990", "-9.99", "-9.990")
  )

  x |>
    tidyr::pivot_longer(cols = all_of(month.abb),
                        names_to = "mon_str", values_to = value_name) |>
    dplyr::mutate(month = match(mon_str, month.abb)) |>
    dplyr::select(-mon_str)
}

oni <- read_psl_12("ONI.txt", "ENSO")
pdo <- read_psl_12("PDO.txt", "PDO")
amo <- read_psl_12("AMO.txt", "AMO")

# --- Read CPC daily index files ---
read_cpc_daily <- function(path, value_name) {
  d <- read_csv(path, comment = "#", show_col_types = FALSE)
  nm <- names(d)
  ycol <- nm[grepl("year", nm, ignore.case = TRUE)][1]
  mcol <- nm[grepl("^mo|month", nm, ignore.case = TRUE)][1]
  dcol <- nm[grepl("^dy|day",   nm, ignore.case = TRUE)][1]
  vcol <- nm[grepl("index|value|nao|ao", nm, ignore.case = TRUE)][1]
  stopifnot(all(!is.na(c(ycol, mcol, dcol, vcol))))

  d |> transmute(
    year  = .data[[ycol]],
    month = .data[[mcol]],
    day   = .data[[dcol]],
    !!value_name := .data[[vcol]]
  )
}

nao <- read_cpc_daily("norm.daily.nao.cdas.z500.19500101_current.csv", "NAO")
ao  <- read_cpc_daily("norm.daily.ao.cdas.z1000.19500101_current.csv",  "AO")

# --- Merge all indices ---
clim_day <- full_join(nao, ao, by = c("year", "month", "day"))
monthly_list <- list(amo, oni, pdo)

clim_idx <- reduce(
  monthly_list,
  .init = clim_day,
  .f = ~ left_join(.x, .y, by = c("year", "month"))
) |> arrange(year, month, day)

# --- Merge with precipitation data ---
load("clim.Rdata")
clim_idx$Date <- as.Date(paste(clim_idx$year, clim_idx$month,
                                clim_idx$day, sep = "-"))

data <- data %>%
  left_join(clim_idx %>% select(Date, NAO, AO, AMO, ENSO, PDO), by = "Date")

save(data, file = "eqte.Rdata")
