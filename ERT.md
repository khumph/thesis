---
title: "ERT"
output: html_document
---

# ERT Algoirthm (Geurts, 2006)
M = number of trees

**`Build_an_extra_tree_ensemble(S).`**  
*Input*: a training set $S$.   
*Output*: a tree ensemble $\mathcal{T} = \{t_1, \ldots, t_M\}.$  

- For $i = 1$ to $M$  
    * Generate a tree: $t_i$ = **`Build_an_extra_tree(S)`**;  
- Return $\mathcal{T}$.


**`Build_an_extra_tree(S)`**.
*Input*: a training set $S$.
*Output*: a tree $t$.
– Return a leaf labeled by class frequencies (or average output, in regression) in $S$ if

1. $|S| < n_min$, or
2. all candidate attributes are constant in $S$, or 
3. the output variable is constant in $S$
    
Otherwise:  

1. Select randomly $K$ attributes, {a1 ,...,aK }, without replacement, among all (non constant in S) candidate attributes;
2.Generate $K$ splits{s1,...,sK},wheresi =Pick a random split(S, ai), ∀i = 1,...,K;
3. Select a split s∗ such that Score(s∗, S) = maxi=1,...,K Score(si , S);
4. Split S into subsets Sl and Sr according to the test s∗;
5. Build $t_l$ = **`Build_an_extra_tree`($S_l$)**. and $t_r$ = **`Build_an_extra_tree(S)`**. from these subsets;
6. Create a node with the split s∗, attach tl and tr as left and right subtrees of this node and return the resulting tree t.

**`Pick_a_random_split(S, a)`**  
*Input*: a training set S and an attribute a.  
*Output*: a split.  
– If the attribute a is numerical

* Compute the maximal and minimal value of a in S, denoted respectively by amS in and amS ax ; * Draw a cut-point ac uniformly in [amS in , amS ax ];
* Return the split [a < ac].

– If the attribute a is categorical (denote by A its set of possible values):

* Compute AS the subset of A of values of a that appear in S;
* Randomly draw a proper non empty subset A1 of AS and a subset A2 of A\AS ; * Return the split [a \in A1 \cup A2].




