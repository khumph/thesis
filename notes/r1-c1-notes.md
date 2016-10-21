# Chapter 1: The reinforcement learning problem

reinforcement learning, is much more focused on goal-directed learning from interaction than are other approaches to machine learning.

like many topics whose names end with “ing,” such as machine learning and mountaineering, is simultaneously a problem, a class of solution methods that work well on the class of problems, and the field that studies these problems and their solution methods.

Reinforcement learning problems:

- involve learning what to do (how to map situations to actions) so as to maximize a numerical reward signal. (much like other machine learning areas)

Unique to RL-
- are closed-loop:
	the learning system’s actions influence its later inputs

- do not involve direct instructions re: what aciton to take: 
	learner is not told which actions to take, but instead must discover which actions yield the most reward by trying them out.

- consequences of actions, including reward signals, play out over extended time periods
	In the most interesting and challenging cases, actions may affect not only the immediate reward but also the next situation and, through that, all subsequent rewards.

basic idea - simply to capture the most important aspects of the real problem facing a learning agent interacting with its environment to achieve a goal.


agent/learner: 
- can sense state of environment
- can take actions to affect the state
- has a goal(s) relating to environment

Any method that is well suited to solving such problems we consider to be a reinforcement learning method.


Reinforcement learning is different from supervised learning,

Supervised learning is learning from a training set of labeled examples (the correct actions the learner should take) provided by a knowledgable external supervisor.

In interactive problems it is often impractical to obtain examples of desired behavior that are both correct and representative of all the situations in which the agent has to act. In uncharted territory—where one would expect learning to be most beneficial—an agent must be able to learn from its own experience.

Also different from unsupervised learning: reinforcement learning is trying to maximize a reward signal instead of trying to find hidden structure.


Unique problem in RL is the trade off between:
Exploitation: use what already known to get reward
Exploration: try new things to see what might work better

Reinforcement learning explicitly considers the whole problem of a goal-directed agent interacting with an uncertain environment.

### 1.2 Examples
Reinforcement learning all involves: 

- interaction between an active decision-making agent and environment
- agent seeks a goal despite uncertainty about environment
- agents actions affect environment, which affect options/opportunities to agent later
- need to take into account later consequences of actions (foresight/planning)

### 1.3 Elements of reinforcement learning
beyond agent & environment
1. policy 
  - way agent behaves at a given time
2. reward signal 
  - immediate feedback (reward) from environment (pleasure or pain) agent's objective is to maximize this reward
  - agent cannot (directly) alter reward signal process
3. value function
  - specifies what's good in long run 
  - value of state = total amount of reward agent can get from that state
  - harder to determine than rewards
4. model of environment (optionally)

model based methods - use models and planning
model-free methods - explicitly trial and error 
- almost opposites


Most important component of most reinforcement learning algorithms is efficiently estimating values

#### evoloutionary methods 
- create many different learners with many different polices, pick the few who create maximum rewards (rather than value functions)

- may offer useful solutions in circumstances esp when learning agent cannot accurately sense environment

- however ignores lots of useful structure:
  - policy searching for is function from states to actions
  - don't notice which states agent passes through
  - or which action agent selects

(note that these may be misleading if agent cannot sense environment well)

#### policy gradient methods
- also don't use value function, but can incorporate - not a sharp distincion
- search in set of parameters, estimate change in parameters that most rapidly improves policy performance
- these estimates are made in real time (not evolutionary time), so take advantage of details of individual interactions

#### temporal difference method

If we let s denote the state before the greedy move, and s' the state after the move, then the update to the estimated value of s, denoted V (s), can be written as
V(s) <- V(s) + alpha (V(s') - V(s)),
where alpha is a small positive fraction called the step-size parameter, which influences the rate of learning. 

called temporal difference because its changes are based on a difference V(s') - V(s) between estimates at two different times

reduce step size over time, or if not reduced to 0, can adapt to slowly changing player

#### contrast evolutionary with value function methods
evo method plays lots with fixed policy (or simulates with model opponent), policy change only after many games, only final outcome used (what happens during game ignored)
value funct methods allow individual states to be evaluated -> takes advantage of info available during game

#### Key features of reinforcement learning methods
- emphasis on learning while interacting with environ (other player)
- clear goal
- correct behavior requires foresight (takes into account delayed effects)
can achieve effects of planning and lookahead without a model of opponent or explicity search over future states/actions


How well a reinforcement learning system can work in problems with such large state sets is intimately tied to how appropriately it can generalize from past experience. 
-> It is in this role that we have the greatest need for supervised learning methods with reinforcement learning.

rl applied even when part of state is hidden or two states appear the same

rl can use models of consequences of actions, but they're not required

some rl methods do not need any model of environ (model-free)
- can be helpful if bottleneck is constructing an accurate enough model of environ
- also important building block for model-based methods

rl can be applied at different "levels" - hierarchical learning systems

# 1.6 Summary
Reinforcement learning is a computational approach to understanding and automating goal-directed learning and decision-making.
different from others
- emphasis on learning by agent from direct interaction with environ w/o examples or complete models of environ

uses formal framework defining interaction between agent and environ: states, actions, rewards.
- simple way of representing essential features of AI problem: cause and effect, uncertainty/nondetermininism and explicit goals
value fucntions important for efficient search in space of policies - distinguishes reinforcement learning from evo methods

# 1.7 History

We define a reinforcement learning method as any ef- fective way of solving reinforcement learning problems

credit assignment problem: how do you ditribute credit for a success among the many decisiion that produced it


missing from supervised learning - drive to achieve some result from the environment, to control the environment toward desired ends and away from undesired ends. This is the essential idea of trial-and-error learning.

local reinforcement - subcomponents ofan overall learning system could reinforce one another

generalized reinforcement - every component (neuron) views all of its inputs in reinforcemtn terms excittatory inputs as rewards, inhibatory inputs ad punishments

uncanny similarity between the behavior of temporal-di↵erence algorithms and the activity of dopamine producing neurons in the brain