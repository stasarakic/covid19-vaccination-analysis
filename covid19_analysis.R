# COVID-19: poredjenje epidemioloskih pokazatelja
# Autorke: Stasa Rakic i Dunja Roza
# Izvor: Our World in Data
# Rekonstrukcija R koda iz dodatka zajednickog izvestaja.

# 1. UCITAVANJE PODATAKA
# Izaberite originalni CSV fajl iz Our World in Data.
podaci <- read.csv(file.choose(), stringsAsFactors = FALSE)

potrebne_kolone <- c(
  "country", "date", "people_fully_vaccinated_per_hundred",
  "new_deaths_smoothed", "new_cases_smoothed", "new_deaths", "new_cases"
)
nedostaju <- setdiff(potrebne_kolone, names(podaci))
if (length(nedostaju) > 0) {
  stop(paste("CSV nema sledece kolone:", paste(nedostaju, collapse = ", ")))
}

podaci$date <- as.Date(podaci$date)

# 2. IZDVAJANJE DRZAVA
drzave <- subset(podaci, country %in% c("Australia", "United Kingdom", "Russia"))
australija <- subset(drzave, country == "Australia")
uk <- subset(drzave, country == "United Kingdom")
rusija <- subset(drzave, country == "Russia")

# 3. FUNKCIJA ZA ANALIZU PO DRZAVI
analiziraj_drzavu <- function(d, naziv) {
  d <- d[order(d$date), ]

  # Prvi datum na kojem je procenat potpuno vakcinisanih >= 50
  dostignut_prag <- !is.na(d$people_fully_vaccinated_per_hundred) &
    d$people_fully_vaccinated_per_hundred >= 50
  if (!any(dostignut_prag)) {
    stop(paste("Prag od 50% nije pronadjen za:", naziv))
  }
  datum50 <- min(d$date[dostignut_prag])

  # Po 180 dana pre i posle datuma dostizanja 50%.
  # Sam datum praga ulazi u period POSLE.
  pre <- d[d$date >= datum50 - 180 & d$date < datum50, ]
  posle <- d[d$date >= datum50 & d$date < datum50 + 180, ]

  # Dnevni odnos sedmodnevno upros. broja umrlih i zarazenih
  odnos_pre <- as.numeric(pre$new_deaths_smoothed) /
    as.numeric(pre$new_cases_smoothed)
  odnos_posle <- as.numeric(posle$new_deaths_smoothed) /
    as.numeric(posle$new_cases_smoothed)
  odnos_pre <- odnos_pre[is.finite(odnos_pre)]
  odnos_posle <- odnos_posle[is.finite(odnos_posle)]

  if (length(odnos_pre) == 0 || length(odnos_posle) == 0) {
    stop(paste("Nema dovoljno validnih dnevnih odnosa za:", naziv))
  }

  # Deskriptivna statistika
  statistike <- function(x) {
    c(Mean = mean(x), Median = median(x), SD = sd(x),
      Min = min(x), Max = max(x))
  }
  deskriptiva <- cbind(
    Pre = statistike(odnos_pre),
    Posle = statistike(odnos_posle)
  )

  # Vilkoksonov test za dva nezavisna uzorka
  vilkokson <- wilcox.test(
    odnos_pre, odnos_posle,
    alternative = "greater", exact = FALSE
  )

  # Zbir dnevnih novih slucajeva i smrti za oba perioda
  umrli_pre <- sum(pre$new_deaths, na.rm = TRUE)
  zarazeni_pre <- sum(pre$new_cases, na.rm = TRUE)
  umrli_posle <- sum(posle$new_deaths, na.rm = TRUE)
  zarazeni_posle <- sum(posle$new_cases, na.rm = TRUE)

  # Proporcije preminulih medju prijavljenim slucajevima
  p_pre <- umrli_pre / zarazeni_pre
  p_posle <- umrli_posle / zarazeni_posle

  # Test jednakosti dve proporcije
  test_proporcija <- prop.test(
    x = c(umrli_pre, umrli_posle),
    n = c(zarazeni_pre, zarazeni_posle),
    alternative = "greater",
    correct = FALSE
  )

  # Prikaz rezultata
  cat("\n====================================\n")
  cat("DRZAVA:", naziv, "\n")
  cat("Datum dostizanja 50%:", as.character(datum50), "\n")
  cat("Pre:", as.character(min(pre$date)), "do",
      as.character(max(pre$date)), "\n")
  cat("Posle:", as.character(min(posle$date)), "do",
      as.character(max(posle$date)), "\n\n")
  cat("DESKRIPTIVNA STATISTIKA\n")
  print(round(deskriptiva, 4))
  cat("\nVILKOKSONOV TEST\n")
  print(vilkokson)
  cat("\nTEST PROPORCIJA\n")
  print(test_proporcija)
  cat("\nProporcija pre:", round(p_pre, 5), "\n")
  cat("Proporcija posle:", round(p_posle, 5), "\n")

  # Grafik 1: boxplot
  boxplot(
    odnos_pre, odnos_posle,
    names = c("Pre 50%", "Posle 50%"),
    main = paste(naziv, "- dnevni odnos umrlih i zarazenih"),
    ylab = "new_deaths_smoothed / new_cases_smoothed",
    col = c("lightblue", "lightgreen")
  )

  # Grafik 2: proporcije
  barplot(
    c(p_pre, p_posle),
    names.arg = c("Pre 50%", "Posle 50%"),
    main = paste(naziv, "- proporcije"),
    ylab = "Odnos ukupnog broja umrlih i zarazenih",
    ylim = c(0, max(p_pre, p_posle) * 1.2),
    col = c("lightblue", "lightgreen")
  )

  invisible(list(
    datum50 = datum50,
    deskriptiva = deskriptiva,
    vilkokson = vilkokson,
    test_proporcija = test_proporcija,
    p_pre = p_pre,
    p_posle = p_posle
  ))
}

# 4. ANALIZA SVE TRI DRZAVE
rezultat_aus <- analiziraj_drzavu(australija, "Australija")
rezultat_uk <- analiziraj_drzavu(uk, "Velika Britanija")
rezultat_rus <- analiziraj_drzavu(rusija, "Rusija")

# Napomena: analiza je opservaciona. Dobijene razlike ne dokazuju
# uzrocni efekat vakcinacije; vremenska zavisnost i drugi faktori
# ogranicavaju tumacenje testova.
