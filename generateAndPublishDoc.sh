#!/bin/bash
pushd .
cd $HOME/svnwork/OPAL/opal-Doc/doc/OPAL/user_guide
make 
makeindex opal_user_guide
make
cp opal_user_guide.pdf /afs/psi.ch/project/amas/www/docs/opal/opal_user_guide-1.1.6.pdf 
popd
exit

doxygen
cp -r doc/html /afs/psi.ch/project/amas/www/docs/opal/

