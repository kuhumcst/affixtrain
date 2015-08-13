Description of workflow that takes a tiny training data set containing words,
lemmas, word classes and some data that is ignored in this example.
The training data set is clustered, so homographs and words that seem to have
the same lemma are gathered in the same cluster. (This clustering is not
essential but improves the quality of the validation step.)

First affixtrain is used to generate binary flexrules in subfolders
    0 1 2 3 4 5
(The folder names indicate the pruning level of the flexrules in that folder).

Then you go to the folder with the unpruned flexrules, where you run affixtrain
to create a pretty printed version of the flexrules and a flexrule file in a
format that Bracmat (https://github.com/BartJongejan/Bracmat) can read.

Go back to the original folder and run Bracmat.

Load a script that reads the bracmat code that utilises the flexrules (in
Bracmat format) to lemmatise input words. The script also reads the original
training data to retrieve the word classes. By lemmatising all training words,
the script now knows which rules are used for which word classes. This 
information is added to the answer when asking to lemmatise a word.

These are the steps:

prompt> affixtrain -i short
prompt> cd 0
prompt> affixtrain -b flexrules.short_XC
prompt> cd ..
prompt> bracmat
{?} get$get$"lemmaVal.bra"
{?} !r

Now enter words and get lemmas and word classes in return.
To stop, enter an empty line to get out of the script's
loop. To exit from Bracmat, enter ) after the prompt.
