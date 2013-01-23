METH=git
METH=https

if [ ! -d letterfunc ]; then
    mkdir letterfunc
    cd letterfunc
    git init
    git remote add origin $METH://github.com/kuhumcst/letterfunc.git
    cd ..
fi
cd letterfunc
git pull origin master
cd ..

if [ ! -d affixtrain ]; then
    mkdir affixtrain
    cd affixtrain
    git init
    git remote add origin $METH://github.com/kuhumcst/affixtrain.git
    cd ..
fi
cd affixtrain
git pull origin master
cd src
make all
cd ..
cd ..
