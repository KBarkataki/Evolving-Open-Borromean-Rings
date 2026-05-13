# Evolving-Open-Borromean-Rings
This package keeps a record of how to reproduce the results presented in Example 4.6 of BarkatakiK,PanagiotouE. 2022TheJonespolynomialofcollectionsofopencurvesin3-space.Proc. R. Soc. A 478: 20220302.https://doi.org/10.1098/rspa.2022.0302

borr_plot.py			:			code to plot a sequence of open borromean rings as they transition continuously from the open configuration to the closed configuration. Also stores the Jones polynomial for each borromean ring in the sequence in a file entitled 'borr_lists.txt'
									The folder 'values_in_paper' contains the sequence of 8 borromen rings used to create Figure 9 in the paper. The polynomial expressions are stored in 'poly_values.txt'

polynomial_plot.py		:	plots a graph of Jones polynomial vs Jones variable for each borromean ring in a chosen sequence.

JONES.py				:			MAIN module for calculating Jones polynomial of collections of open and/or closed curves in 3-space.

functions.py			:			module of auxillary functions used within JONES.py

# Whenever you access this code package 
python3 -m venv borr_env

source borr_env/bin/activate

pip install -r requirements.txt

python borr_plot.py
