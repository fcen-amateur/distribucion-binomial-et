library(tidyverse)

# Funciones de soporte
cumverosimilitud <- function(X, tita) {
  dbinom(X, size = tita$k, prob = tita$p) %>% log %>% cumsum
}

verosimilitud <- function(X, tita) {
  dbinom(X, size = tita$k, prob = tita$p) %>% log %>% sum
}

muhat <- function(X) {
  mean(X)
}

sigma2hat <- function(X) {
  mean((X - mean(X))^2)
}

criterio_levin <- function(X) {
  # Si este cociente es <= 1, la verosimilitud tiene máximo en tita(+Inf, 0)
  muhat(X) / sigma2hat(X)
}

# Estimadores para k, p y tita, según EMM y EMV
fun_estimadores <- list(
  k = list(
    mom1 = function(X, tita) {
      khat <- mean(X) / tita$p
      return(tita(round(khat), NA))
    },
    mom2 = function(X, tita) {
      p <- tita$p
      khat <- (sqrt(p ^ 2 - 2 * p + 4 * mean(X ^ 2) + 1) + p - 1) / (2 * p)
      return(tita(round(khat), NA))
    },
    maxv = mv_k <- function(X, tita) {
      max_vero <- -Inf
      k_actual <- max(X)
      while (verosimilitud(X, tita(k_actual, tita$p)) > max_vero) {
        max_vero <- verosimilitud(X, tita(k_actual, tita$p))
        k_actual <- k_actual + 1
      }
      return(tita(k_actual, NA))
    }
  ),
  p = list(
    mom1 = function(X, tita) {
      phat <- mean(X) / tita$k
      return(tita(NA, phat))
    },
    mom2 = function(X, tita) {
      k <- tita$k
      phat <- (k + sqrt(k ^ 2 + 4 * k * (k - 1) * mean(X ^ 2))) / (2 * k * (k - 1))
      return(tita(NA, phat))
    },
    maxv = function(X, tita) {
      phat <- mean(X) / tita$k
      return(tita(NA, phat))
    }
  ),
  tita = list(
    mom = function(X, tita) {
      k <- round(mean(X) ^ 2 / (mean(X) ^ 2 + mean(X) - mean(X ^ 2)))
      p <- mean(X) / k
      return(tita(k, p))
    },
    maxv = function(X, tita) {
      k_infinito <- criterio_levin(X) <= 1
      if (k_infinito) {
        return(tita(+Inf, 0))
      }
      media_X <- mean(X)
      tita_actual <- tita(k = max(X), p = media_X / max(X))
      max_vero <- -Inf
      while(verosimilitud(X, tita_actual) > max_vero) {
        max_vero <- verosimilitud(X, tita_actual)
        tita_actual$k <- tita_actual$k + 1
        tita_actual$p <- media_X / tita_actual$k
      }
      return(tita_actual)
    }
  )
)

estimadores <- tribble(
  ~param, ~metodo,
  "k",    "mom1",
  "k",    "mom2",
  "k",    "maxv",
  "p",    "mom1",
  "p",    "mom2",
  "p",    "maxv",
  "tita",    "mom",
  "tita",    "maxv"
)

estimar <- function(param, metodo, X, tita) {
  fun_estimadores[[param]][[metodo]](X, tita)
}

samplear <- function(n, tita) {
  rbinom(n, size = tita$k, prob = tita$p)
}

tita <- function(k = NA, p = NA) {
  list(k = k, p = p)
}

computar_estimadores <- function(estimadores, X, tita) {
  estimadores %>%
    mutate(
      titahat = map2(param, metodo, estimar, X = X, tita = tita))
}

ronda <- function(estimadores, rango_n, rango_k, rango_p) {
  
  n <- round(runif(1, rango_n[1], rango_n[2]))
  k <- round(runif(1, rango_k[1], rango_k[2]))
  p <- runif(1, rango_p[1], rango_p[2])
  
  tita0 <- tita(k, p)
  X <- samplear(n, tita0)
  
  computar_estimadores(estimadores, X, tita0) %>%
    transmute(param, metodo,
             khat = map_dbl(titahat, "k"),
             phat = map_dbl(titahat, "p"),
             k = k,
             p = p,
             n = n
    )
}

principal <- function(cant_sims, estimadores, rango_n, rango_k, rango_p,
                      ruta_output = "simulaciones") {
  contador <- 1
  wd <- getwd()
  while (contador < cant_sims) {
    resultados <- ronda(estimadores, rango_k, rango_n, rango_p)
    write_csv(resultados,
              path = paste0(wd, "/", ruta_output, "/", 
                            digest::digest(resultados), ".csv"))
    contador <- contador + 1
  }
}

rango_n <- c(1, 300)
rango_k <- c(1, 1000)
rango_p <- c(0, 1)
cant_sims <- 100

# Muestra tomada del paper de Olkin, Petkau & Zidek
tita_olkin <- tita(75, 0.32)
X_olkin <- c(16, 18, 22, 25, 27)


#principal(cant_sims, estimadores, rango_n, rango_k, rango_p, ruta)

