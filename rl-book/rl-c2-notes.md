banit problem - special case of rl problem wher there is only a single state

rl uses training infor that evaluates the actions taken rather than instructs by giving correct actions

evaluative feedback - indicates how good the action taken is (but not whether it is the best/worst action)
  - basis for function optimization (including evo methods)

instructive feedback - indicates the correct action to take, independent of the action actually taken. 
  - basis of supervised learning
  
nonassociative setting - does not involve learning to act in more than one situation

# 2.1 k-armed bandit problem
- faced repeatedly with a choice among k options (actions)
- after each choice, you receive a numerical reward from a stationary prob dist depending on the action
- objective: maximize reward over some time period (e.g. 100 action selections called time steps)

name comes from a nickname for a slot machine (one armed bandit) but in this case you have k choices. A doctor treating successive seriously ill patients with experimental treatments is ssame prob.

value of action - expected reward given that action is selected

greedy actions = the action at any time step with highest estimated value

If you take a greedy action, you are exploiting your current knowledge of the values of the actions

If you don't take the greedy action, you are exploring (enables improvement of estimate of that action)

one often refers to "conflict" between exploration and exploitation (b/c can't do both in single action)

finding optimal balance is difficult in practice (requires strong assumptions)

we'll consider simple cases

# 2.2 action-value methods

sample-average method - natural way to estimate true value of action is to add up rewards and divide by the number of times action was taken (sample mean) 

greedy selection method - always pick action with highest estimated value

epsilon-greedy methods - pick action with highest expected value with probability 1 - epsilon but select at random (uniformly) from all actions w.p. epsilon.

10-armed testbed

nonstationary - the true values of the actions changed over time

effective nonstationarity is the case most commonly encountered in reinforcement learning

Even if the underlying task is stationary and deterministic, the learner faces a set of banditlike decision tasks each of which changes over time due to the learning process itself. 

# incremental implementation


NewEstimate <- OldEstimate + StepSize[Target - OldEstimate]

[Target - OldEstimate] = error in estimate, reduced by takng a step towards the target (nth reward in the above)

step-size parameter denoted by alpha 
or alpha_t(a) [step taking action a at time t]

# tracking a nonstationary problem
weight recent events more heavily than long-past ones

popular method: use a constant step size parameter alpha.

creates a weighted average of rewards, more recent weighted heavier

this method sometimes called an exponential, recency-weighted average.

alpha_n(a) = 1/n results in simple-average (converges to true action values by LLN)

if inf sum of alpha_n(a)s is inf and inf sum of squared alpha_n(a)s is finite, then results converge to true action values

fixed alpha doesn't meet this (which is good for notstationary)

# 2.5 optimistic initial values

above methods depend on initial aciton value estimates (Q_1(a)). Called biased by their initial estiates.

bias disappears for sample-average methods once all acitons have been taken, but not for constant alpha (though dec over time)

in practice this bias not a problem, but initial estimates become parameters that must be picked by user, but can be an easy way to supply prior knowledge.

optimistic initial values: 
can also use to encourage exploration: set values high
useful on stationary problems, but not well suited to nonstationary (drive for exploration is temporary)

in general, initial state is of less interest in nonstationary

# 2.6 Upper-Confidence-Bound Action Selection

it would be better to explore nongreedy options by the chance of them being optimal rather than at random

eqn 2.8

performs well but is ore difficult to extend to general settings (nonstationary problems)

# 2.8 gradient bandits
numerical preference H_t(a) for each action a. Larger preference, more often action is taken. only relative preference matters no interpretatio re: reward.

action probs defined by eqn 2.9

gradient ascent - each preference H_t(a) would be incrementing proportional to the increments effect on performance (2.11)
since we don't know  true value (but 2.10 has same expected value, called stochasitc gradient ascent - has nice convergence properties)

# 2.9 associative search (contextual bandits)

several k-armed bandit tasks that change randomly as play - would appear as a single nonstationary k-armed bandit whose true action values change randomly, but you have some clue re: which bandit you're facing

associatve serach - involves both trial-and-error (search for best actions) and association of these actions with the situations in which they are best (somtimes called contextual bandits)

intermediate between k-armed bandit and general rl problem

They are like the full reinforcement learning problem in that they involve learning a policy, but like our version of the k-armed bandit problem in that each action a↵ects only the immediate reward.

# 2.9 Summary


