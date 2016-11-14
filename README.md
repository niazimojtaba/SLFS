# SLFS
Author : Mojtaba Niazi (If you want to use this package please email to me mojtaba.niazi.095@gmail.com)
This code refer to SLFS algorithm. You can have source article from "http://www.sciencedirect.com/science/article/pii/S0952197616301105". 
<b>Abstract</b>

Feature selection is an important task in many problems occurring in pattern recognition, bioinformatics, machine learning and data mining applications. The feature selection approach enables us to reduce the computation burden and the falling accuracy effect of dealing with huge number of features in typical learning problems. There is a variety of techniques for feature selection in supervised learning problems based on different selection metrics. In this paper, we propose a novel unified framework for feature selection built on the graphical models and information theoretic tools. The proposed approach exploits the structure learning among features to select more relevant and less redundant features to the predictive modeling problem according to a primary novel likelihood based criterion. In line with the selection of the optimal subset of features through the proposed method, it provides us the Bayesian network classifier without the additional cost of model training on the selected subset of features. The optimal properties of our method are established through empirical studies and computational complexity analysis. Furthermore the proposed approach is evaluated on a bunch of benchmark datasets based on the well-known classification algorithms. Extensive experiments confirm the significant improvement of the proposed approach compared to the earlier works.

<b>Keywords</b>

    Feature selection; Supervised learning; Relevant features; Mutual information; Structure learning; Graphical models

How to run SLFS algorithm on linux based os

-First you should compile "main.cpp" file. Write below command in your terminal in main.cpp directory.

<code>g++ main.cpp -o SLFS

-What things you should notice about input file?
SLFS uses 9-fold cross validation for calculating accuray, So you should combine all data in a text file as input. In the input file each rows of data are samples from your data, So each sample contions features. The last feature is class vaiable. Finally data should be seperated with space.

-How to run?

Write below command in your execution folder:

<code>SLFS inputFile.txt outputFile.txt</code>

-The outputFile contions results from each fold, The first line is accuray value, the third line consist of features indexes, 
the forth line shows mutual information value with class, and the fifth line is number of selected features with class variable. these lines repeated for each part of cross validation.

-How can you have the Bayesian Classifier?

You should uncomment a part of code in bfs_ function that write adjacent matrix on adj.txt file.

 
