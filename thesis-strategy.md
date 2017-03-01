#### goals: generalize CRT to more covariates, more models
- [ ] add noise covariates (don't interact with treatment), method should remove them
- [ ] add variables that interact with treatment, method should include them: be able to identify subgroups that respond better
    - [ ] generalize to simple cutoffs interacting with treatment
- [ ] reward depending on previous stage
- [ ] apply to real data with transitions

# possible tasks
Add some additional interacting variables
actually make plots based off predicted versus actual for last time point?
second pass Q learning and A learning paper

# possibilities
## split into two projects
### 1. an investigation into high dimensional data in a generic smart
- several treatments, subgroups, re-randomize nonresponders
- just one continuous outcome (e.g. weight lost)
- easier to make misspecified models?
_I feel like knowing what to do with high dimensional data or how to identify subgroups is less interesting in our setting because it's not a design anyone actually uses_

### 2. Picking dose based on survival outcome
- make setup more realistic (see below)
- show that it works across different parameterizations
then incorporate high dimensions, subgroups

## variable importance
- [x] [_in principle done through varImp, but not yet widely implemented_]
- [ ] identify specifically what variables are important for subgroups _do we even have subgroups in the original formulation?_

### add another continuous drug with a different mechanism?

### change randomization to be closer to a SMART
- [ ] More efficient way to do it besides random doses? - doesn't seem very realistic that people would consent to a trial carried out the way in the paper
- [ ] Do we run into problems that the assumption of the DTRs we're considering may not actually be represented in the data, i.e, all possible trajectories aren't represented?
- [ ] make the randomization more explicitly analogous to a SMART? seems like you lose some of the nice features of a SMART design in the Clinical reinforcement trial, e.g. non-responders get a new treatment, etc

### Compare to other strategies actually used in practice
- [ ] Compare to play the winner type strategies?
- [ ] Actually test the "widely used" approach of max dose, then time off, then max dose?

## models
_potential models that will do variable selection "automatically" and allow nonlinear relationships_

#### models that do both
- [x] GCV MARS
- [x] ERT
- [x] rpart
- [ ] rf via ranger (can force dose to be included)

#### models that do variable selection but not nonlinearities
- [x] LASSO

#### models that do nonlinearities but not variable selection
- [ ] SVR _does this do variable selection?_
- [ ] GAM (_the restricted cubic spline model fails for high dimensionality df > n_)

#### models that do neither
- [ ] OLS

#### tune models to get better performance
_currently tuning parameters are picked via default caret tuning: 25 bootstrap samples on 3 different combinations of parameters (for each Q-function)_

## survival outcome as reward
_probably also useful for VADT_
- [x] log of survival time as reward - no censoring
- [x] **[Just use survival time as reward]** - **evaluate sensitivity of proximal to distal rewards** - if proximal rewards are not entirely predictive of distal rewards: seems like same sorts of issues with surrogate outcomes can arise with picking utilities, method probably really depends on intermediate outcomes being prognostic of real outcome of interest - look at misspecification of proximal rewards? e.g. toxicity doesn't really matter at all in mortality

- the original way rewards defined contradicts the ideas in S&B (don't reward intermediate outcomes)

refs: Q-learning with censored data paper, follow up NSCLC paper, DTR book (references one other paper)

## More scenarios

#### different noise scenarios
- [x] 10 variables have an effect out of 100
- [ ] 20 variables have an effect out of 100
- [ ] 50
- [ ] different correlations between noise for above?

#### different interaction/subgroup scenarios

#### different subgroups and noise scenarios

##### time-varying noise variables

#### misspecify model
some strange function that interacts with treatment e.g. polynomial try to show that it is robust (neural network type as in elem stat learn MARS section)

## modify transition functions to avoid perfect balance on survival

- [ ] go back and figure out where doses are the same, and choose the doses that end up with best (even if tied)?

_there is a perfect balance except for when either W or M is 0:_
<iframe src="https://www.desmos.com/calculator/esw7qu2rx5?embed" width="500px" height="500px" style="border: 1px solid #ccc" frameBorder=0 />

Attempt at modification:
<iframe src="https://www.desmos.com/calculator/fq2m2fdqls?embed" width="500px" height="500px" style="border: 1px solid #ccc" frameBorder=0 />

- [ ] Change the params on mus on survival? generate from uniform distribution?

- [x] [_just remove balance, make actual max for survival_] make more sophisticated way to choose max for ties in "best" based off of sum of tumor mass and toxicity?

**make a new project that should be more or less identical to paper**?
- [ ] replicate ERT results in CRT paper
  - [ ] how is my ERT different from their's?
- [ ] replicate SVR results in CRT paper
  - [ ] get SVR working with one stage of thesis data

## testing/refactoring
- [ ] extract repeated sequences as functions

### unit tests for each piece

# other ideas
- [ ] See how sample size effects results?
- [ ] average over many simulations? _Would be more accurate, but not very realistic, as in reality you only observe one trial_
- [ ] average test set over (e.g. 10) replicates? _Not sure this is necessary since the expected survival time is used in test set. But would give different starting conditions for the 200 patients_


***

# future directions
- [ ] Sample size formulae for CRT with survival outcomes
- [ ] confidence intervals for CRT with survival outcomes
- [ ] incorporating adaptive elements into a SMART? - "in some settings the incorporation of adaptive elements into a SMART design is possible (Thall et al. 2002, Thall & Wathen 2005), how to achieve optimal incorporation is an open question that warrants further research."

***

# Qs to answer
- [ ] How do these Bayesian Thall strategies ignore treatment heterogeneity? Do the Berry ones too?
- [ ] How is the Q learning and and A learning methods for dynamic treatment regimes relevant to response depending on previous stage?
- [ ] How does the RLT design help with drug discovery?
- [ ] search literature for sequential design of experiments?
- [ ] How to have a single terminal utility/reward instead of sum of stages in Q learning? - DTR review says this can also be done
- [ ] What does assuming sufficient regularity mean? - DTR review paper
- [ ] What are nonregular asymptotics?
