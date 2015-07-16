affixtrain
==========

Using supervised learning, create a set of affix rules for use by the CSTlemma
lemmatiser.

Training takes place in three stages:

  * In the first stage the program tries to find a set of six parameters that
    optimize the tree produced for the current language. In practice, the best
    way of optimization is finding a (local) minimum for the number of tree
    nodes. Other optimizations can give marginally better results.
  * In the second stage the program trains trees for increasing numbers of
    training examples, starting with a very small percentage and ending with
    almost all available training examples. The remaining examples are used for
    testing. Depending on the size of the training or test set, the procedure
    is repeated for different random samples. For each training set size the
    precision is estimated, together with standard deviation. The parameters
    are computed for the power law that gives the expected number of rules for
    a given number of training words. The lower the exponent, the better is the
    generalization effect. For each training set size the rules are also pruned
    a number of times. For each pruning level the precision is computed anew.
    Often, the best rules (fewest wrongly lemmatised words) are obtained by
    pruning the rules that are supported by fewer than three examples.
    A tabular report is added (as a comment) to a parameter file that can be
    used for future reference and re-training.
  * In the third stage all available training words are used to produce
    production-ready rules, again in several pruning levels. It is up to the
    user to decide the pruning level to use. Unpruned rules lemmatise all
    training words correctly, but may give more erroneous lemmas for out-of-
    vocabulary words. Pruned rules may give better results for OOV words, but
    will not lemmatize all training words correctly. To cope with this
    disadvantage you should also provide a binary dictionary to the lemmatiser.
    Use cstlemma -h to see how to do that.
    
Empty lines in the training/testing data are interpreted as cluster separators.
If the data has no empty lines between non-empty lines, the training and
testing occurs on a line-by-line basis, but if there are such empty lines,
training and testing occurs on a cluster-by-cluster basis. For example, by
collecting homographs in clusters and defining all non-ambiguous full forms as
one-line clusters, testing with 'OOV' words (that is, words that were not used
during training) will result in more realistic estimates of how well the rules
are able to spot and lemmatize ambiguous full forms.
    
Notice that the whole process easily can take many days, even a couple of
weeks, to run.
    
This version still contains a lot of "dead wood" and confusing naming. We work
on that.

Bart Jongejan, April 21, 2015