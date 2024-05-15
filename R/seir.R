#' SEIR
#'
#' This function runs an SEIR model simulation and plots the results.
#'
#' @return Returns a plot of the SEIR model simulation.
#' @export
#' @importFrom ggplot2 ggplot geom_line aes labs scale_color_manual theme element_text ggsave
#' @importFrom deSolve ode
#' @importFrom stats time
#' @importFrom utils tail

model_and_plot <- function() {
  # 定义模型函数
  model <- function(t, y, param) {
    S <- y[1]
    E <- y[2]
    I1 <- y[3]
    I2 <- y[4]
    R <- y[5]
    N <- param["N"]
    beta <- param["beta"]
    mu <- param["mu"]
    gamma1 <- param["gamma1"]
    gamma2 <- param["gamma2"]
    lamda <- param["lamda"]
    sigma <- param["sigma"] #I2感染率

    dSt <- mu * (N - S) - beta * S * I1/N
    dEt <- beta * S * I1/N - mu * E - lamda * E
    dI1t <- lamda * E - (mu + gamma1 + sigma) * I1
    dI2t <- sigma * I1 - (mu + gamma2) * I2
    dRt <- gamma1 * I1 + gamma2 * I2 - mu * R
    outcome <- c(dSt, dEt, dI1t, dI2t, dRt)
    list(outcome)
  }

  # 参数
  times <- seq(0, 156, by = 1/7)
  param <- c(mu = 0.000, lamda = 0.03, beta = 4, gamma1 = 0.1, gamma2 = 0.05, sigma = 0.01, N = 1)
  init <- c(S = 0.9999, E = 0.00008, I1 = 0.00001, I2 = 0.00001, R = 0)

  result <- deSolve::ode(y = init, times = times, func = model, parms = param)
  result <- as.data.frame(result)
  names(result) <- c("time", "S", "E", "I1", "I2", "R")
  tail(round(result, 6), 10)
  utils::globalVariables(c("S", "E", "I1", "I2", "R"))

  # 绘图
  seriplot <- ggplot2::ggplot(data = result, ggplot2::aes(x = time)) +
    ggplot2::geom_line(ggplot2::aes(y = S, color = "S"), lwd = 2) +
    ggplot2::geom_line(ggplot2::aes(y = E, color = "E"), lwd = 2) +
    ggplot2::geom_line(ggplot2::aes(y = I1, color = "I1"), lwd = 2) +
    ggplot2::geom_line(ggplot2::aes(y = I2, color = "I2"), lwd = 2) +
    ggplot2::geom_line(ggplot2::aes(y = R, color = "R"), lwd = 2) +
    ggplot2::labs(x = "Time", y = "Ratio", title = "SEIR Model Dynamics") +
    ggplot2::scale_color_manual(name = "SEIRD",
                       values = c(S = "blue", E = "purple", I1 = "red", I2 = "orange", R = "green")) +
    ggplot2::theme(legend.title = ggplot2::element_text(size = 12),
          legend.text = ggplot2::element_text(size = 10))
  print(seriplot)
  # Save the plot to files
  ggplot2::ggsave(filename = "seir.pdf", plot = seriplot, width = 7, height = 6)
  ggplot2::ggsave(filename = "seir.svg", plot = seriplot, width = 7, height = 6)
}

# 调用函数，运行模型并绘图
model_and_plot()
