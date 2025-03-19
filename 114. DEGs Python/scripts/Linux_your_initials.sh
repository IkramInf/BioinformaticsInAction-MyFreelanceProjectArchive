#!/bin/bash

# concatenating the gene expression matrices for variety X 
cat VarX/*.csv | sort -u > matrixX.csv
sed 's/-a/Rep.1/g; s/-b/Rep.2/g; s/-c/Rep.3/g' VarX_Header.csv > HeaderX.csv
cat HeaderX.csv matrixX.csv > VarX_Matrix.csv

# concatenating the gene expression matrices for variety Y
cat VarY/*.csv | sort -u > matrixY.csv
sed 's/-a/Rep.1/g; s/-b/Rep.2/g; s/-c/Rep.3/g' VarY_Header.csv > HeaderY.csv
cat HeaderY.csv matrixY.csv > VarY_Matrix.csv
