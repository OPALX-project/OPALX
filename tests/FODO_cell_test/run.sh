set -e

python calc.py

cd opal
opal test.in --info 10
cd ..

python test.py
python hd5.py