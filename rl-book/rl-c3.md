# 3.1 The Agent–Environment Interface

The reinforcement learning problem is meant to be a straightforward framing of the problem of learning from interaction to achieve a goal.

agent - The learner and decision-maker

environment - The thing the agent interacts with, comprising everything outside the agent

A complete specification of an environment defines a task, one instance of the reinforcement learning problem.

agent and environment interact at a sequence of discrete time steps

at each step time t,
  1. agent receives some representation of the environment's state ($S_t \in \mathcal{S}$)
  2. on basis of state, agent selects an action ($A_t \in \mathcal{A}(S_t)$)
  3. one time step later, the agent receives a numerical reward ($R_{t+1} \in \mathcal{R} \subset \mathbb{R}$)
  4. then the agent finds itself in a new state, $S_{t+1}$



# 3.2 Goals and Rewards

reward hypothesis:

> "That all of what we mean by goals and purposes can be well thought of as the maximization of the expected value of the cumulative sum of a received scalar signal (called reward)."

Reward signal - used to tell agent what we want achieved, not how to achieve it.

agent's goal: maximize cumulative reward

tell agent what we want:
critical that maximizing rewards actually achieves our aim; rewards truly indicate what is desired

Not how to achieve it:

  - don't reward actions taken on the way to achieving goal
  - this could lead to maximizing subgoals while not achieving overall objective (e.g. chess)

things that may seem internal to agent can be considered part of environment (e.g. limbs, power sources)

convenient to place limit of agent at the limit of its direct control (not necessarily its body)
  - agent's ultimate goal should be under imperfect control
    - can't simply decree it gets a reward (b/c rewards come from environment)

# 3.3 Returns

In general we seek to maximize expected return (G_t) - some function of rewards, e.g. the sum of rewards at each time step

## Episodic tasks
- where there is a natural notion of a final time step, followed by a reset

episode - a natural repeated subsequence (e.g. plays of a game, trips through a maze)

each episode ends in a terminal state, then resets

Set of all nonterminal states denoted $\matcal{S}$

set of all states plus terminal state, $\mathcal{S}^+$

## Continuing tasks

interaction does not have defined episodes, but continues indefinitely

Since plain sum of rewards could now be infinite, seeks to maximize discounted sum:

\[
G_t \doteq R_{t+1} + \gamma R_{t+2} + \gamma^2 R_{t+3}  + \cdots = \sum_{i = 0}^{\infty} \gamma^{i} R_{t + i + 1}
\]

Where $\gamma$ is a parameter 0 \leq \gamma \leq 1 called the discount rate
