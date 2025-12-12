## Adaptive System Identification Using LS, RLS, and Kalman Filter

This report presents an extensive simulation study on **adaptive system identification** of a coupled mass–spring–damper system with translational and rotational dynamics.

The system is modeled in continuous time, discretized using an appropriate sampling strategy, and identified using both **offline and online estimation techniques** under a wide range of operating conditions.
<img width="855" height="178" alt="image" src="https://github.com/user-attachments/assets/cba3cf47-e636-4f9a-8e56-173beb094175" />

### Topics Covered
- Continuous-time and discrete-time state-space modeling
- Transfer function derivation and discretization
- Least Squares (LS) system identification
- Recursive Least Squares (RLS) and Extended LS (ELS)
- Forgetting-factor RLS (λ-RLS)
- Covariance resetting strategies for time-varying parameters
- Kalman Filter (KF) based parameter estimation
- Identification under feedback and instability
- Nonlinear system identification

### Excitation and Noise Analysis
The estimation performance is evaluated using multiple excitation signals:
- Impulse, step, ramp
- Sinusoidal and multi-sinusoidal inputs
- White noise and colored noise

The effect of **persistent excitation (PE) order** on parameter convergence is explicitly analyzed. Results show strong dependence of estimation accuracy on excitation richness and noise color.

### Model Order Study
The report includes a systematic comparison between:
- Exact-order models
- Under-parameterized models
- Over-parameterized models

The impact of model mismatch on estimation quality and output tracking is demonstrated through numerical simulations.

### Key Observations
- LS performs well under sufficient PE and white noise but degrades under colored noise
- RLS-based methods improve robustness for online and time-varying systems
- Forgetting factors and covariance resetting significantly enhance tracking of parameter changes
- Kalman Filter provides competitive performance, particularly in feedback configurations

All results are supported by time-domain responses, parameter convergence plots, and quantitative comparisons.
