# terms
reinforcement learning - lookahead notion

sparsity - a small set of variables completely convey the true signal

improved performance over random forests in high dimensional settings

asymptotic properties - properties as sample size goes to infinity

greedy algorithm - making the locally optimal choice at each stage

nearest neighbor

checkerboard - little or no marginal effects, but strong joint effects

# abstract
implements reinforcement learning at each selection of a splitting variable 

splits on the variable that brings the greatest future improvement in later splits, rather than the one with largest marginal effect from the immediate split, -> uses the available samples in a more efficient way.

approach enables linear combination cuts at little extra computational cost

variable muting procedure - eliminates noise variables - also uses reinforcement learning
 - eliminates from future consideration: so toward terminal nodes when sample size is small, noise variables aren't considered

# 1. intro

#### drawbacks to rf 
- large p small n - random feature selection in rf - most important “random” component of random forests and the driving force be- hind the improvement from a bagging predictor, not well suited to sparse data

- little opportunity to consider a strong variable as the splitting rule and would also lead to bias in the variable importance measures

- a large number of predictors causes overfitting toward terminal nodes where the sample size is small and prevents the effect of strong variables from being fully explored.



#### Biau 2012
- addressed concerns re: rf with a special type of prf where strong vars have higher prob of selection
  - model forces all or most splits on strong variables
- if this prob can be properly chosen, convergence rate depends only on number of strong variables

However,

- prob of using strong var to split at internal node depends on within-node data (or independent set of within-node samples). 
  - samp size small near terminal nodes unlikely to behave well for entire tree: in sim by Biau, sample size < 25 prob of using strong var can be very low
- marginal comparisons of splitting vars (esp in high dim settings) can potentially fail to identify strong vars (some may have little marginal but strong joint effects)


#### general features:
- choose variables for each split that will bring largest return from future branching splits (rather than immediate consquences via marginal effects)
  - break any hidden structure and avoid inconsistency by forcing splits on strong variables even if they do not show any marginal effect

- progressively mute noise variables as go deeper down tree
  - even as the sample size decreases rapidly toward a terminal node, the strong variable(s) can still be properly identified from the reduced space

- enables linear combination splitting rules at very little extra computational cost.
  - "The linear combination split constructed from the strong variables gains efficiency when there is a local linear structure and helps preserve randomness under this somewhat greedy approach to splitting variable selection."
  
#### Consequences
- forces the splits to concentrate only on the p1 strong variables at the early stage of the tree construction while also reducing the number of candidate variables gradually toward terminal nodes. 
  - more sparse tree structure
    - consider less variables esp toward terminal nodes
- convergrnce under certian assumptions depends on the number of strong predictors not total number of predictors

# 2. proposed method
In short, the proposed RLT model is a traditional random forests model with a special type of splitting variable selection and noise variable muting.

#### Checkerboard example 
- 2 strong vars, only joinly predictive of outcome
- trad rf won't pick 2 strongs over others
- if we know in advance that splitting on either of strongs will be advantageous, we could force split on them

How to identify strongs?

- fit a pilot model evaluate the potential contribution of each variable.
(which is embedded at each internal node, and thus will be called an embedded model throughout)

- then split node using identified most important variable(s)

When doing this recursively for each daughter node, we can focus the splits on the variables that will very likely lead to a tree yielding the smallest prediction error in the long run.

problems with greedy splitting var selection:
ample size shrinks as we move toward a terminal node, it becomes increasingly difficult to identify the important variables regardless of what embedded model we are using

extreme concentration on the strong variables could lead to highly correlated trees even when bootstrapping is employed.

solutions proposed: 
variable muting procedure to counter the first drawback and use linear combination splits to introduce extra randomness. 


embedded model: a model fitted to internal node A data, if A is root, no variables will be muted, but if A isn't then some will be (more as you move down tree).
 
use modification of ET (extremely randomized trees) fitting each tree with bootstrapped sample  (ETs are attractive because they can be faster - reduce computational cost), but other methods (PRF, or rf) can be used

variable importance  uses breiman (rf) definition  

#### 2.6 muting 
muted set - set we've muted p_d

protected set - vars that aren't muted. Anything split on is in protected set

minimal number p0 of vars beyond which no more are muted. updated after split is done

@root:
pick vars with p0 highest importance at root, those are protected
pick vars with pd lowest importance, those go in muted

@internal node
pick var(s) split on (argmax VI) (more than 1 if lincom is used), add lowest pd VI to muted set (of those that arent' in protected)

NB: pd important tuning parameter, can vary based on how many vars are not in muted set

