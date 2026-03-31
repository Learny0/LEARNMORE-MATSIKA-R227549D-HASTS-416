# ============================================================
# HASTS 416 – Tutorial 1 in R 
# LEARNMORE MATSIKA R227549D
# ============================================================

#install.packages(c("markovchain", "igraph", "expm"))

library(markovchain)
library(igraph)
library(expm)

options(repr.plot.width = 10, repr.plot.height = 6)

# ============================================================
# HELPER FUNCTION
# ============================================================

state_period <- function(P, s, states, max_n = 100) {
  i <- which(states == s)
  visits <- which(sapply(1:max_n, function(n) (P %^% n)[i, i] > 0))
  if (length(visits) == 0) return(Inf)
  Reduce(function(a, b) { while (b != 0) { tmp <- b; b <- a %% b; a <- tmp }; a },
         diff(c(0, visits)))
}

# ============================================================
# QUESTION A1 - 5-State Markov Chain
# ============================================================

P1 <- matrix(c(
  1,0,0,0,0,
  0.5,0,0,0,0.5,
  0.2,0,0,0,0.8,
  0,0,1,0,0,
  0,0,0,1,0),
  nrow=5, byrow=TRUE)

states1 <- c("S1","S2","S3","S4","S5")
rownames(P1) <- colnames(P1) <- states1
mc1 <- new("markovchain", transitionMatrix=P1, states=states1)

# ── A1(a): Diagram & Classification ────────────────────────

par(mar=c(1,1,3,1))
plot(mc1,
     main="A1(a): 5-State Markov Chain",
     layout=layout_in_circle,
     vertex.size=45,
     vertex.label.cex=1.2,
     edge.label.cex=1,
     edge.arrow.size=0.5,
     vertex.color="dodgerblue4")

cat("\n=== A1(a) Classification ===\n")
cat("\nCommunicating classes:\n"); print(communicatingClasses(mc1))
cat("\nRecurrent classes:\n"); print(recurrentClasses(mc1))
cat("\nTransient classes:\n"); print(transientClasses(mc1))
cat("\nAbsorbing states:\n"); print(absorbingStates(mc1))

cat("\nPeriods:\n")
for (s in states1) {
  d <- state_period(P1, s, states1)
  print(paste(s, d))
}

# ── A1(b): Trajectories ─────────────────────────────────────

set.seed(42)
n_steps <- 30
state_num <- function(s) match(s, states1)

traj <- replicate(3, {
  s0 <- sample(states1,1)
  c(s0, rmarkovchain(n_steps, mc1, t0=s0))
}, simplify=FALSE)

par(mar=c(5,5,4,2))

plot(0:n_steps, sapply(traj[[1]], state_num),
     type="o", pch=16, lwd=2,
     col="blue", ylim=c(1,5),
     yaxt="n",
     xlab="Time Step", ylab="State",
     main="A1(b): Simulated Trajectories")

axis(2, at=1:5, labels=states1, las=1)

lines(0:n_steps, sapply(traj[[2]], state_num),
      type="o", col="red", pch=16, lwd=2)

lines(0:n_steps, sapply(traj[[3]], state_num),
      type="o", col="darkgreen", pch=16, lwd=2)

legend("topright",
       legend=c("Trajectory 1","Trajectory 2","Trajectory 3"),
       col=c("blue","red","darkgreen"),
       lwd=2, pch=16, bty="n")

cat("\nComment: Chains eventually get absorbed into S1.\n")

# ══════════════════════════════════════════════════════════════
# A1(c): STEADY-STATE PROBABILITIES
# ══════════════════════════════════════════════════════════════

cat("\n=== A1(c) Steady-State Probabilities ===\n")

ss1 <- steadyStates(mc1)
cat("\nSteady-state distribution for A1:\n")
print(ss1)

# Check irreducibility
is_irred1 <- is.irreducible(mc1)
cat("\nIs mc1 irreducible?", is_irred1, "\n")

# Check period using period() function
period1 <- period(mc1)
cat("Period of mc1:", period1, "\n")

# For a chain with absorbing states, steady-state is degenerate
if (nrow(ss1) > 0) {
  cat("\nInterpretation:\n")
  cat("The steady-state distribution is [1, 0, 0, 0, 0] - the chain is absorbed in S1.\n")
  cat("State S1 is an absorbing state with probability 1 in the long run.\n")
  cat("All other states (S2, S3, S4, S5) are transient.\n")
  cat("\nIs the chain ergodic?\n")
  cat("ERGDIC = FALSE\n")
  cat("Reason: Chain has an absorbing state (S1), so it is NOT irreducible.\n")
  cat("Ergodic chains must be irreducible, aperiodic, and positive recurrent.\n")
} else {
  cat("\nNo unique steady-state exists.\n")
}

# ── A1(d): Convergence ─────────────────────────────────────

n_time <- 50
init <- rep(1/5,5)
prob_mat <- matrix(0, n_time+1, 5)
prob_mat[1,] <- init

for(t in 1:n_time){
  prob_mat[t+1,] <- prob_mat[t,] %*% P1
}

par(mar=c(5,5,4,5))
matplot(0:n_time, prob_mat,
        type="l", lwd=3, lty=1,
        col=rainbow(5),
        xlab="Time", ylab="Probability",
        main="A1(d): Probability Convergence",
        ylim=c(0,1))

legend("right", legend=states1,
       col=rainbow(5), lwd=3,
       inset=c(-0.2,0), xpd=TRUE, bty="n")

grid()

cat("\nComment on convergence: Probabilities converge rapidly to steady-state.\n")
cat("State S1 probability approaches 1 as n increases.\n")

# ============================================================
# QUESTION A2 - 7-State Markov Chain 
# ============================================================

# ── Define the 7-state transition matrix (from question paper) ──

P2 <- matrix(c(
  # Col1  Col2  Col3  Col4  Col5  Col6  Col7
  0,     1,    0,    0,    0,    0,    0,    # Row 1: S1
  1,     0,    0,    0,    0,    0,    0,    # Row 2: S2
  0,     0,    0,    0.4,  0.2,  0.2,  0.2,  # Row 3: S3
  0,     0,    0,    0,    0.2,  0.4,  0.4,  # Row 4: S4
  0.3,   0,    0,    0.1,  0.3,  0.1,  0.2,  # Row 5: S5
  0,     0,    0,    0.2,  0.2,  0.3,  0.3,  # Row 6: S6
  0,     0,    0,    0.5,  0.2,  0.2,  0.1   # Row 7: S7
), nrow=7, byrow=TRUE)

states2 <- c("S1","S2","S3","S4","S5","S6","S7")
rownames(P2) <- colnames(P2) <- states2
mc2 <- new("markovchain", transitionMatrix=P2, states=states2)

# ── A2(a): Diagram ──────────────────────────────────────────

par(mar = c(2, 2, 3, 2))
plot(mc2,
     main               = "A2(a): 7-State Markov Chain",
     edge.arrow.size    = 0.35,
     edge.label.cex     = 0.7,
     edge.curved        = 0.35,
     vertex.color       = "darkorange",
     vertex.label.color = "white",
     vertex.label.cex   = 0.95,
     vertex.size        = 30,
     layout             = layout_in_circle)

# ── A2(b): Classification ──────────────────────────────────

cat("\n=== A2(b) Classification ===\n")

cat("\nCommunicating classes:\n"); print(communicatingClasses(mc2))
cat("\nRecurrent classes:\n"); print(recurrentClasses(mc2))
cat("\nTransient classes:\n"); print(transientClasses(mc2))
cat("\nAbsorbing states:\n"); print(absorbingStates(mc2))
cat("\nAbsorbing (in usual sense):\n")
abs_states <- which(apply(P2, 1, function(row) sum(row) == 1 && row[which(row == 1)] == 1))
cat("States that return to themselves with probability 1:",
    if(length(abs_states) > 0) states2[abs_states] else "None", "\n")

cat("\nReflecting states (self-loop with probability > 0):\n")
ref_states <- which(diag(P2) > 0)
cat("States with self-loops:",
    if(length(ref_states) > 0) states2[ref_states] else "None", "\n")

cat("\nPeriods:\n")
for (s in states2) {
  d <- state_period(P2, s, states2)
  print(paste(s, ": period =", d))
}

# ── A2(c): Trajectories ────────────────────────────────────

set.seed(2024)
n_steps2 <- 50
state_num2 <- function(s) match(s, states2)

traj2 <- replicate(2, {
  s0 <- sample(states2,1)
  c(s0, rmarkovchain(n_steps2, mc2, t0=s0))
}, simplify=FALSE)

par(mar=c(5,5,4,2))
plot(0:n_steps2, sapply(traj2[[1]], state_num2),
     type="o", col="blue", pch=16,
     ylim=c(1,7), yaxt="n",
     xlab="Time Step", ylab="State",
     main="A2(c): Two Simulated Trajectories")

axis(2, at=1:7, labels=states2, las=1)

lines(0:n_steps2, sapply(traj2[[2]], state_num2),
      type="o", col="red", pch=16)

legend("topright", legend=c("Trajectory 1","Trajectory 2"),
       col=c("blue","red"), lwd=2, pch=16)

cat("\nA2(c) Comment: Trajectories starting from different states tend to\n")
cat("move towards the communicating classes. States S1 and S2 form a\n")
cat("2-state periodic class (period 2). States S3-S7 form another class.\n")

# ── A2(d): Limiting Distribution ───────────────────────────

cat("\n=== A2(d) Limiting Probabilities ===\n")

is_irred2 <- is.irreducible(mc2)
cat("\nIs mc2 irreducible?", is_irred2, "\n")

# Get period for each state
period2 <- period(mc2)
cat("Period of mc2:", period2, "\n")

ss2 <- steadyStates(mc2)
cat("\nLimiting/Stationary probabilities:\n")
print(ss2)

if (nrow(ss2) > 0) {
  cat("\nLimiting distribution exists.\n")
} else {
  cat("\nNo unique limiting distribution - multiple recurrent classes.\n")
}

cat("\nErgodic check:\n")
cat("For a chain to be ergodic, it must be:\n")
cat("  1. Irreducible: all states communicate\n")
cat("  2. Aperiodic: period = 1\n")
cat("  3. Positive recurrent: finite expected return time\n")

if (!is_irred2) {
  cat("\nCONCLUSION: The chain is NOT ergodic.\n")
  cat("Reason: Chain is NOT irreducible (has multiple communicating classes).\n")
  cat("Communicating classes: {S1,S2} and {S3,S4,S5,S6,S7}\n")
} else if (period2 > 1) {
  cat("\nCONCLUSION: The chain is NOT ergodic.\n")
  cat("Reason: Chain is periodic with period", period2, "\n")
} else {
  cat("\nCONCLUSION: The chain is ergodic.\n")
}

# ============================================================
# QUESTION A3 - Traffic Markov Chain
# ============================================================

P_day <- matrix(c(0.4,0.4,0.2,
                  0.3,0.4,0.3,
                  0,0.1,0.9),
                nrow=3, byrow=TRUE)

P_peak <- matrix(c(0.1,0.5,0.4,
                   0.1,0.3,0.6,
                   0,0.1,0.9),
                 nrow=3, byrow=TRUE)

mat_power <- function(M,n){
  result <- diag(nrow(M))
  for(i in 1:n) result <- result %*% M
  result
}

# ── A3(a): Distribution at 6PM ─────────────────────────────

pi0 <- c(1,0,0)  # Start with Light state at 1PM

# 1PM-4PM: 3 hours = 180 min / 20 min intervals = 9 steps
# 4PM-6PM: 2 hours = 120 min / 20 min intervals = 6 steps

pi_6pm <- pi0 %*% mat_power(P_day,9) %*% mat_power(P_peak,6)

cat("\n=== A3(a) Distribution at 6PM ===\n")
cat("\nStarting from Light state at 1PM, distribution at 6PM:\n")
print(pi_6pm)
colnames(pi_6pm) <- c("Light","Heavy","Jammed")
print(pi_6pm)

cat("\nInterpretation:\n")
cat("After 5 hours (15 time intervals), starting from Light:\n")
cat("  - Probability of Light:", round(pi_6pm[1,1], 4), "\n")
cat("  - Probability of Heavy:", round(pi_6pm[1,2], 4), "\n")
cat("  - Probability of Jammed:", round(pi_6pm[1,3], 4), "\n")

# ── A3(b): Simulation Verification ─────────────────────────

set.seed(123)
N <- 10000

simulate_trajectory <- function(){
  state <- 1  # Start at Light
  # 1PM-4PM (9 steps)
  for(i in 1:9) state <- sample(1:3,1,prob=P_day[state,])
  # 4PM-6PM (6 steps)
  for(i in 1:6) state <- sample(1:3,1,prob=P_peak[state,])
  state
}

res <- replicate(N, simulate_trajectory())
emp <- table(res)/N

# Convert to proper format
emp_vec <- c(emp[1], emp[2], emp[3])
names(emp_vec) <- c("Light","Heavy","Jammed")

comparison <- rbind(Analytical=pi_6pm[1,], Simulation=emp_vec)
cat("\n=== A3(b) Simulation Comparison (N=10000) ===\n")
print(comparison)

cat("\nSimulation verifies analytical calculation.\n")
cat("Differences are due to sampling error.\n")

# ── Plot Comparison ────────────────────────────────────────

par(mar=c(5,5,4,2))
bp <- barplot(comparison,
              beside=TRUE,
              col=c("blue","red"),
              ylim=c(0,1),
              main="A3: Analytical vs Simulation Distribution at 6PM",
              ylab="Probability",
              xlab="Traffic State")

legend("topleft",
       legend=c("Analytical","Simulation (N=10000)"),
       fill=c("blue","red"),
       bty="n")

text(bp, comparison + 0.02,
     labels=round(comparison,3),
     cex=0.9)

# ============================================================
# END OF COMPLETE SOLUTIONS
# ============================================================
