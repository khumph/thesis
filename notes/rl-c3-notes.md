---
output:
  html_document:
    self_contained: yes
---


# 3.1 The Agent–Environment Interface

The reinforcement learning problem is meant to be a straightforward framing of the problem of learning from interaction to achieve a goal.

agent - The learner and decision-maker

environment - The thing the agent interacts with, comprising everything outside the agent

A complete specification of an environment defines a task, one instance of the reinforcement learning problem.

agent and environment interact at a sequence of discrete time steps

at each step time t,

1. agent receives some representation of the environment's state ($S_t \in \mathcal{S}$)
2. on basis of state, agent selects an action $(A_t \in \mathcal{A}(S_t))$
3. one time step later, the agent receives a numerical reward ($R_{t+1} \in \mathcal{R} \subset \mathbb{R}$)
4. then the agent finds itself in a new state, $S_{t+1}$



# 3.2 Goals and Rewards

Reward hypothesis:

> "That all of what we mean by goals and purposes can be well thought of as the maximization of the expected value of the cumulative sum of a received scalar signal (called reward)."

Reward signal - used to tell agent what we want achieved, **not how to achieve it**.

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

Set of all nonterminal states denoted $\mathcal{S}$

set of all states plus terminal state, $\mathcal{S}^+$

## Continuing tasks

interaction does not have defined episodes, but continues indefinitely

Since plain sum of rewards could now be infinite, seeks to maximize discounted sum:

\[
G_{t} \doteq R_{t+1} + \gamma R_{t+2} + \gamma^2 R_{t+3}  + \cdots = \sum_{i = 0}^{\infty} \gamma^{i} R_{t + i + 1}
\]

Where $\gamma$ is a parameter $0 \leq \gamma \leq 1$ called the discount rate


# 3.4 Unified Notation for Episodic and Continuing Tasks
Use $S_t$, etc for both episodic and continuing tasks

We almost always refer to a single episode, or state that something is true for all episodes when considering episodic tasks, so instead of referring to $S_{t,i}$ (the state at time $t$ for episode # $i$), we drop the $i$.

consider the terminal state to be an absorbing state

absorbing state: state that transitions only to itself and generates only rewards of zero


# 3.5 Markov Property


# 3.6 Markov Decision Processes

Most theory implicitly assumes the environment is a finite MDP


expected rewards for state–action pairs

state-transition probabilities

expected rewards for state–action–next-state triples,


All you need to know to

Given any state and action s and a, the probability of each possible pair of next state and reward, s',r is completely specify the dynamics of a finite MDP


# 3.7 Value functions

Almost all RL algorithms involve estimating value functions

Value function:
functions of states (or state-action pairs) that estimate how good it is for the agent to be in a given state (or how good it is to perform a given action in a given state)

"how good" is defined in terms of future rewards that can be expected: expected return

Since rewards depend on what actions a agent takes, value functions are defined with respect to particular polices

#### What is a policy $(\pi)$ in reinforcement learning?
A policy, $\pi$, is a mapping from each state $s \in \mathcal{S}$, and action $a \in \mathcal{A}(s)$, to the probability, $\pi(a \mid s)$, of taking action a when in state $s$.

The policy, $\pi$, tells us the probability of taking action $a$ when in state $s$.

#### What is the value of a state $s$ under a policy $\pi$?
Informally, the value of state $s$ under policy $\pi$ (denoted $v_{\pi}(s)$) is the expected return when starting in $s$ and following $\pi$ thereafter.

#### What is the value function for MDPs?

\[
v_{\pi}(s) \doteq \mathbb{E}_{\pi}[G_{t} \mid S_{t} = s] =
  \mathbb{E}_{\pi}\left[ \sum_{k=0}^{\infty} \gamma^{k} R_{t+k+1} \mid S_{t} = s \right]
\]

called the *state-value function for policy $\pi$*


#### What is the value of an $a$ under policy $\pi$?
Denoted $q_{\pi}(s,a)$
The expected return starting from $s$, taking action $a$, and then following $\pi$ thereafter


\[
q_{\pi}(s,a) \doteq \mathbb{E}_{\pi}[G_{t} \mid S_{t} = s, A_{t} = a] =
  \mathbb{E}_{\pi}\left[
    \sum_{k=0}^{\infty} \gamma^{k} R_{t+k+1} \mid S_{t} = s, A_{t} = a
  \right]
\]

called the *action-value function for policy $\pi$*

#### What are Monte Carlo methods?

averaging over many random samples of actual results

if an agent follows policy $\pi$ and maintains an average, for each state encountered, of the actual returns that have followed that state, then the average will converge to the state’s value, v_{\pi}(s), as the number of times that state is encountered approaches infinity. If separate averages are kept for each action taken in a state, then these averages will similarly converge to the action values, q$\pi$ (s, a).

if there are very many states, then it may not be practical to keep separate averages for each state individually.

Instead, the agent would have to maintain v$\pi$ and q$\pi$ as parameterized functions and adjust the parameters to better match the observed returns. This can also produce accurate estimates, although much depends on the nature of the parameterized function approximator (Chapter 9).


A fundamental property of value functions used throughout reinforcement learning and dynamic programming is that they satisfy particular recursive relationships.

What recursive relationships do value functions satisfy?

bellman eqn p 17-18

backup diagrams

# 3.8 Optimal Value Functions

These operations transfer value information back to a state (or a state– action pair) from its successor states (or state–action pairs).

"Value functions define a partial ordering over policies" ordering over??

For finite MDPs, A policy $\pi$ is better than another policy $\pi$' if its expected return is greater than or equal to $\pi$' for all states:

$\pi $\pi$' if and only if v_{\pi}(s) \geq v{\pi'}(s) for all s \in \mathcal{S}.

There is always one policy that is better than or equal to all others
