dynamic treatment regime - list of decision rules that formalizes the treatment decisions made for each patient based on their evolving situation

Q denotes "quality" in Q-learning

1. Category: What type of paper is this? A description and comparison of research methods
2. Context: Lots of other papers involving Q and A learning. Which theoretical bases were used to analyze the problem?
3. Correctness: Do the assumptions appear to be valid? Sure?
4. Contributions: What are the paper’s main contributions? Reference, resource for determining if Q or A learning should be used in a given scenario
5. Clarity: Is the paper well written? yes.


k = 1, ..., K decision points
Y = outcome (large values preferred)
Ω = superpopulation of patients
ω ∈ Ω = patient from this population
S = covariates
A_k = treatment options
̄a_k = (a_1, ..., a_k) possible treatment history
W^* = potential outcomes including Y\*(\bar{a}_k), and covariate information between decisions, S\*


What assumptions are required to use observed data to estimate an optimal regime?
