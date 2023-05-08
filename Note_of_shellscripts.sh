# Note of shell scripts

## virtualenv of python3 is installed with pip.
## virtualenv is different from venv
#python -m pip install --user virtualenv
virtualenv medaka --python=python3 --prompt "(medaka)"
##activate the virtial env
. medaka/bin/activate
pip install --upgrade pip
pip install medaka

## deactivate the virtialenv
deactivate