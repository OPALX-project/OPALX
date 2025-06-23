set -e

python calc.py

cd opal
opal test.in
cd ..

python test.py