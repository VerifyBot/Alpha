# <img style="position:relative;top:10px" src="https://i.imgur.com/aZ1ZtyJ.png" float="right" width="50" height="50"> Alpha

Nir Yona -- Alpha bioinformatics project using **Python** @ https://www.cs.tau.ac.il/~tamirtul


*<p>Full project paper **SOON!**</p>*

# 🪁 Tools
* Python (pandas, sklearn, matplotlib/plotly, numpy, pickle, etc.)
* Jupyter Notebook & Google Colab
* ViennaRNA, Chimera UGEM

# 🎯 Goal

This project focused on creating features and traing a model that will predict re-initiation chances for E-Coli operons.

## Methods

### Data

The data was provided from a library created [in this paper](https://pubmed.ncbi.nlm.nih.gov/32973167/) - `~13,000`
synthetic variants. The project focused on the randomic part of the data, as shown in this image:

<img src="https://i.imgur.com/uFpfFbt.png" width="300">


### Generating Features

To predict the re-initiation chances best, we created textual features that describe the sequence and the structure of
the sequence (like occurances of different nucs in certain places and reading layers, etc) and
biological features known to be related to protein expression and therefore to re-initiation chances (
like [CAI](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-9-65), [tAI](https://academic.oup.com/bioinformatics/article/33/4/589/2593585), [Folding levels](https://pubmed.ncbi.nlm.nih.gov/32973167/)
and more)
After the features creation a Multiple Regression model starts training.
At the end we try the feature on E.Coli and test p-value to find out the meaningful features.

### Training

We trained a Multiple Regression model on a library variants with their re-initiation chances.
The model was trained on 60% of the data, tested on 20% and validated on the last 20%.
After the training and testing we tried to optimize the features, and by that get the features that
affect the re-initiation chances the most. We use cross validation to make sure we don't overfit the data.

### Analyzing the Features on E.Coli Operons

To make sure the chosen features are really significant, we tried them on endogenous E.Coli operons.
We created premutations of the operons, and finally calculated p-value for each feature against the library operons.

## 🎯 Results

### Multiple Regression Model

The model was trained over and over adding more feautures every run, and the prediction results were as follows:

<img src="https://i.imgur.com/0FUgQxU.png" width="500">

### Features Analysis

For the significant features we set the p-value border to 0.05 as a significant feature requirement.

<img src="https://i.imgur.com/ILU3pw7.png" height="250">

<br>

<br>

**The results conclude that the amount of appearances of the `TT` dinucleotide affects the re-initiation chances in a
linear way the most out of all the other features calculated!**


The results are backed up by other papers showing similar relationships (like [Yang et al., 2019](https://academic.oup.com/nar/article/47/17/9243/5549713) or [Nie, Wu, & Zhang, 2006](https://pubmed.ncbi.nlm.nih.gov/17028312/))

## 🐉 Run

Create a virtual environment and install the requirements from `requirements.txt`,
then run the desired notebook with jupyter notebook or the script with python.

```bash
py -m pip install requirements.txt
jupyter notebook
```
