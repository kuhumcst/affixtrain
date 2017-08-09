Description of workflow that takes a tiny training data set containing words,
lemmas, word classes and some data that is ignored in this example.
The training data set is clustered, so homographs and words that seem to have
the same lemma are gathered in the same cluster. (This clustering is not
essential but improves the quality of the validation step.)

First affixtrain is used to generate flexrules in a binary format in the
subfolders
    0 1 2 3 4 5
The folder names indicate the pruning level of the flexrules in that folder.

To start with, train flexrules while ignoring the word class information in
the training data, so only using the full forms in the first column and the
lemmas in the secons column. For this, specify the option -n FLOO (or just
-n FL)

If you have a good POS-tagger for your language, you may want to train
flexrules for each word class or tag. Use the -n parameter to specify in
which column full form, lemma and word class (or tag) are to be found. But
even if you train flexrules for each word class or tag, you should also train
flexrules while ignoring word class and tag. CSTlemma uses these neutral
flexrules for words that have a word class or tag that has no representation in
the training data.

If you want to get some insight in the created rules, you can proceed as 
follows. 

Go to the folder with the unpruned flexrules, where you run affixtrain
to create a pretty printed version of the flexrules and a flexrule file in a
format that Bracmat (https://github.com/BartJongejan/Bracmat) can read.

Then go back to the original folder and run Bracmat.

Load a script that reads the bracmat code that utilises the flexrules (in
Bracmat format) to lemmatise input words. The script also reads the original
training data to retrieve the word classes. By lemmatising all training words,
the script now knows which rules are used for which word classes and how often.
This information is added to the answer when asking to lemmatise a word.

These are the steps:

prompt> affixtrain -i short -n FLOO
prompt> cd 0
prompt> affixtrain -b flexrules.short_XC
prompt> cd ..
prompt> bracmat
{?} get$"lemmaVal.bra"
{?} !r

Now enter words (for exampe eighties and reducing) and get lemmas and a list
of possible word classes in return.
To stop, enter an empty line to get out of the script's
loop. To exit from Bracmat, enter ) after the prompt.
