# scDHA_ac
The scDHA_ac can perform cell clustering through unsupervised learning, which is developed by R language. We suggest that you could run it within R-studio.


# Basic Usage
Make sure that R language development environment is ready on your machine.

Clone the repository: `git clone https://github.com/Huyanmei123/AC.git`

The main interface of the function is in adaptive.R

# To run the Goolam example
- New a file, and add the source  `source('adaptive.R')`

- Get the dataset `data<-read.csv("/your path/goolam.csv",header=TRUE)`

- Get its true label  `label<-read.csv("/your_path/goolam_label.csv",header=True)`,`label<-label$cell_type1`

- Log transform the data: `data <- log2(data + 1)`

- Generating clustering result: `result <- adaptive(data, seed = 1)`

- The clustering result can be found here: `cluster <- result$cluster`

- Calculating adjusted Rand Index using aricode package: `aricode::ARI(cluster, label)`

  

