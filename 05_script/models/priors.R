#==============================================================================#
#                                                                              #
#                                Check priors                                  #
#                                                                              #
#==============================================================================#


# U(0, 1) ----------------------------------------------------------------------
q0cubs <- rep(0, 100000)
q1cub <- runif(n = 100000, 0, 1)
q2_3cubs <- runif(n = 100000, 0, 1)

ggplot(data.frame(q1cub), aes(x = q1cub)) +
  geom_histogram(fill = "white", color = "black") +
  xlim(-0.1, 1.1) 

# back-transform
p1cub <- exp(q1cub)/(exp(q0cubs) + exp(q1cub) + exp(q2_3cubs))
p2_3cubs <- exp(q2_3cubs)/(exp(q0cubs) + exp(q1cub) + exp(q2_3cubs))

# Plot
ggplot(data.frame(p1cub), aes(x = p1cub)) +
  geom_histogram(fill = "white", color = "black") +
  xlim(-0.1, 1.1) +
  labs(title = "U(0, 1)")

ggsave("10_meetings/2021-06-XX Meeting with Sarah/unif(0, 1).png", 
       width = 3, height = 3)

# U(0, 10) ----------------------------------------------------------------------
q0cubs <- rep(0, 100000)
q1cub <- runif(n = 100000, 0, 1)
q2_3cubs <- runif(n = 100000, 0, 1)


# back-transform
p1cub <- exp(q1cub)/(exp(q0cubs) + exp(q1cub) + exp(q2_3cubs))
p2_3cubs <- exp(q2_3cubs)/(exp(q0cubs) + exp(q1cub) + exp(q2_3cubs))


# Plot
ggplot(data.frame(p1cub), aes(x = p1cub)) +
  geom_histogram(fill = "white", color = "black") +
  xlim(-0.1, 1.1) +
  labs(title = "U(0, 10)")

ggsave("10_meetings/2021-06-XX Meeting with Sarah/unif(0, 10).png", 
       width = 3, height = 3)


# Norm(0, 10) ------------------------------------------------------------------
q0cubs <- rep(0, 100000)
q1cub <- rnorm(n = 100000, mean = 0, sd = 10)
q2_3cubs <- rnorm(n = 100000, mean = 0, sd = 10)

ggplot(data.frame(q1cub), aes(x = q1cub)) +
  geom_histogram(fill = "white", color = "black") +
  xlim(-0.1, 1.1)

# back-transform
p1cub <- exp(q1cub)/(exp(q0cubs) + exp(q1cub) + exp(q2_3cubs))
p2_3cubs <- exp(q2_3cubs)/(exp(q0cubs) + exp(q1cub) + exp(q2_3cubs))


# Plot
ggplot(data.frame(p1cub), aes(x = p1cub)) +
  geom_histogram(fill = "white", color = "black") +
  xlim(-0.1, 1.1)  +
  labs(title = "N(0, 10)")

ggsave("10_meetings/2021-06-XX Meeting with Sarah/norm(0, 10).png", 
       width = 3, height = 3)


# Norm(0, sd = 1.5) ------------------------------------------------------------------
q0cubs <- rep(0, 100000)
q1cub <- rnorm(n = 100000, mean = 0, sd = 1.5)
q2_3cubs <- rnorm(n = 100000, mean = 0, sd = 1.5)

ggplot(data.frame(q1cub), aes(x = q1cub)) +
  geom_histogram(fill = "white", color = "black") +
  xlim(-0.1, 1.1) 

# back-transform
p1cub <- exp(q1cub)/(exp(q0cubs) + exp(q1cub) + exp(q2_3cubs))
p2_3cubs <- exp(q2_3cubs)/(exp(q0cubs) + exp(q1cub) + exp(q2_3cubs))


# Plot
ggplot(data.frame(p1cub), aes(x = p1cub)) +
  geom_histogram(fill = "white", color = "black") +
  xlim(-0.1, 1.1) +
  labs(title = "N(0, 1.5)") 

ggsave("10_meetings/2021-06-XX Meeting with Sarah/norm(0, 1.5).png", 
       width = 3, height = 3)

# Norm(0, sd = 1.5^2) ------------------------------------------------------------------
q0cubs <- rep(0, 100000)
q1cub <- rnorm(n = 100000, mean = 0, sd = 1.5^2)
q2_3cubs <- rnorm(n = 100000, mean = 0, sd = 1.5^2)

# back-transform
p1cub <- exp(q1cub)/(exp(q0cubs) + exp(q1cub) + exp(q2_3cubs))
p2_3cubs <- exp(q2_3cubs)/(exp(q0cubs) + exp(q1cub) + exp(q2_3cubs))


# Plot
ggplot(data.frame(p1cub), aes(x = p1cub)) +
  geom_histogram(fill = "white", color = "black") +
  xlim(-0.1, 1.1) +
  labs(title = "N(0, 1.5^2)") 

ggsave("10_meetings/2021-06-XX Meeting with Sarah/norm(0, 1.5^2).png", 
       width = 3, height = 3)


# Beta(1, 1) (Same as Unif(1, 1)) ----------------------------------------------
q0cubs <- rep(0, 100000)
q1cub <- rbeta(n = 100000, 1, 1)
q2_3cubs <- rbeta(n = 100000, 1, 1)

# back-transform
p1cub <- exp(q1cub)/(exp(q0cubs) + exp(q1cub) + exp(q2_3cubs))
p2_3cubs <- exp(q2_3cubs)/(exp(q0cubs) + exp(q1cub) + exp(q2_3cubs))


# Plot
ggplot(data.frame(p1cub), aes(x = p1cub)) +
  geom_histogram(fill = "white", color = "black") +
  xlim(-0.1, 1.1) +
  labs(title = "Beta(1, 1)") 

ggsave("10_meetings/2021-06-XX Meeting with Sarah/beta(1, 1).png", 
       width = 3, height = 3)




