

# HASTS 416 - Tutorial 1 Solutions
#R227549D LEARNMORE MATSIKA
#HDSC
---

## **QUESTION A1: 5-State Markov Chain**

### **A1(a) Classification**

**Transition Matrix P:**
```
    S1   S2   S3   S4   S5
S1  1.0  0.0  0.0  0.0  0.0
S2  0.5  0.0  0.0  0.0  0.5
S3  0.2  0.0  0.0  0.0  0.8
S4  0.0  0.0  1.0  0.0  0.0
S5  0.0  0.0  0.0  1.0  0.0
```

**Communicating Classes:**
- `{S1}` - singleton class
- `{S2}` - singleton class
- `{S3, S4, S5}` - communicating class

**Recurrent Classes:** `{S1}` (S1 is the only recurrent class)

**Transient Classes:** `{S2}`, `{S3, S4, S5}`

**Absorbing State:** `S1` (once entered, cannot leave)

**Reflecting States:** None (no self-loops)

**Periods:**
| State | Period |
|-------|--------|
| S1 | 1 (absorbing) |
| S2 | ∞ (transient) |
| S3 | 3 |
| S4 | 3 |
| S5 | 3 |

---

### **A1(b) Trajectories**

Three simulated trajectories starting from random states show the chain eventually getting **absorbed into S1**. All paths flow toward the absorbing state.

---

### **A1(c) Steady-State Probabilities**

**Steady-State Distribution:**
```
S1    S2    S3    S4    S5
1.0   0.0   0.0   0.0   0.0
```

**Is the chain irreducible?** `FALSE`

**Is the chain ergodic?** `FALSE`

**Interpretation:** 
- The chain has an **absorbing state (S1)** with probability 1 in the long run
- All other states are **transient** (probability 0 in steady-state)
- The chain is **NOT ergodic** because it has an absorbing state (not irreducible)

---

### **A1(d) Convergence**

The probability convergence plot shows:
- S1 probability rapidly approaches 1.0
- All other state probabilities converge to 0
- Convergence is fast due to the structure of the chain

---

## **QUESTION A2: 7-State Markov Chain**

### **A2(a) & A2(b) Classification**

**Communicating Classes:**
- `{S1, S2}` - 2-state periodic class
- `{S3}` - singleton transient
- `{S4, S5, S6, S7}` - aperiodic class

**Recurrent Classes:** `{S1, S2}`, `{S4, S5, S6, S7}`

**Transient Classes:** `{S3}`

**Absorbing States:** None

**Reflecting States:** S5, S6, S7 (have self-loops)

**Periods:**
| State | Period |
|-------|--------|
| S1 | 2 |
| S2 | 2 |
| S3 | ∞ (transient) |
| S4 | 1 |
| S5 | 1 |
| S6 | 1 |
| S7 | 1 |

---

### **A2(c) Trajectories**

Trajectories show the chain moving within communicating classes. States S1 and S2 form a **periodic class (period 2)** alternating between them. States S4-S7 form an **aperiodic class**.

---

### **A2(d) Limiting Probabilities**

**Is mc2 irreducible?** `FALSE`

**Stationary Probabilities:**
```
S1    S2    S3    S4    S5    S6    S7
0.5   0.5   0.0   0.0   0.0   0.0   0.0
```

**Is the chain ergodic?** `FALSE`

**Reason:** Chain is NOT irreducible (has multiple communicating classes):
- `{S1, S2}` 
- `{S4, S5, S6, S7}`

---

## **QUESTION A3: Traffic Markov Chain**

### **A3(a) Distribution at 6PM**

Starting from **Light** state at 1PM, distribution at 6PM:

| State | Probability |
|-------|------------|
| **Light** | 0.0148 (1.48%) |
| **Heavy** | 0.1326 (13.26%) |
| **Jammed** | 0.8526 (85.26%) |

**Calculation:** π(6PM) = π₀ × P_day⁹ × P_peak⁶

---

### **A3(b) Simulation Verification (N=10,000)**

| Method | Light | Heavy | Jammed |
|--------|-------|-------|--------|
| Analytical | 0.0148 | 0.1326 | 0.853 |
| Simulation | 0.0146 | 0.1287 | 0.857 |

**Conclusion:** Simulation closely matches analytical results. Differences are due to sampling error (expected with Monte Carlo methods).
