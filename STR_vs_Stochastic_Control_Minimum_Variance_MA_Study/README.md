## Comparison of STR and Stochastic Controllers Under Noise

This report presents a comprehensive simulation-based comparison between
**Self-Tuning Regulators (STR)** and **stochastic control strategies** for discrete-time
dynamic systems under noisy conditions.

The study focuses on **minimum variance (MV)** and **moving average (MA)** controllers
and evaluates their performance in both **adaptive and non-adaptive** settings, considering
direct and indirect controller formulations.

<img width="460" height="328" alt="image" src="https://github.com/user-attachments/assets/b6414aa8-605e-4b5d-b55d-c18afb164f43" />


### Topics Covered
- Dynamic system modeling and transfer function derivation
- Discretization and stability analysis
- Direct and indirect STR with pole–zero cancellation
- STR performance with and without white measurement noise
- Comparison of direct vs indirect STR under noise
- Non-adaptive minimum variance control (colored noise)
- Adaptive minimum variance control using parameter estimation
- Non-adaptive and adaptive moving average control
- Effect of noise modeling on control performance
- Controller behavior under minimum-phase and non-minimum-phase systems

### Noise and Estimation Analysis
The report explicitly considers **colored noise**, incorporating noise dynamics into
controller design. Parameter estimation is performed using **Extended Least Squares (ELS)**,
allowing simultaneous estimation of plant and noise parameters.

The role of noise in enriching the identification signal and accelerating parameter
convergence is analyzed in detail.

### Direct vs Indirect Control Structures
Both direct and indirect formulations are examined across all controllers.
Results show that indirect methods consistently achieve:
- Lower output variance
- Reduced control effort
- Improved robustness under noise

Direct methods exhibit faster initial response but suffer from higher oscillations and
significantly increased control effort, especially in adaptive configurations.

### Minimum Variance vs Moving Average Control
- Minimum variance control performs best for **minimum-phase systems**, achieving the
lowest output variance.
- Moving average control is more suitable for **non-minimum-phase systems**, where zero
elimination is not feasible.
- The trade-off between transient performance and noise suppression is demonstrated
through extensive simulations and statistical comparisons.

### Non-Minimum Phase Systems
A non-minimum-phase system is explicitly constructed to enable stable and meaningful
comparisons. Results show that:
- Minimum variance control struggles in non-minimum-phase scenarios
- Moving average control provides improved stability and variance reduction
- Adaptive MA controllers outperform their non-adaptive counterparts in such cases

### Key Observations
- Indirect controllers outperform direct controllers across all scenarios
- Adaptive controllers improve parameter convergence but increase initial transients
- Noise modeling is essential for effective variance reduction
- MV control is best suited for regulation problems
- MA control is preferable for non-minimum-phase systems
- STR pole placement offers superior transient shaping but weaker noise suppression
