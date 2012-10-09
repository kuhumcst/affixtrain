METH=https

if [ ! -d letterfunc ]; then
    mkdir letterfunc
    cd letterfunc
    git init
    cd ..
fi
cd letterfunc
git remote add origin $METH://github.com/kuhumcst/letterfunc.git
git pull origin master
cd ..

if [ ! -d affixtrain ]; then
    mkdir affixtrain
    cd affixtrain
    git init
    cd ..
fi
cd affixtrain
git remote add origin $METH://github.com/kuhumcst/affixtrain.git
git pull origin master
cd src
make all
cd ..
cd ..

