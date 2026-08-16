# Bayesian Spatio-temporal Regression Model

Simulation code for a Bayesian spatio-temporal model of **group (pooled) testing**
data, using a Gaussian Predictive Process (GPP) approximation with Pólya-Gamma
data augmentation. The study compares three testing protocols for estimating
disease prevalence over space and time.

## Testing scenarios

- **Individual testing (IT)** — each specimen tested separately (baseline).
- **Dorfman testing (DT)** — pools are tested first; specimens in positive pools
  are then retested individually.
- **Master-pool testing (MPT)** — pools tested once, with no individual retesting.
- **Realistic master-pool testing (RMPT)** — pools tested once, with realistic pool sizes, no individual retesting.
