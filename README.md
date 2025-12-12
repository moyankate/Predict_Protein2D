# Protein Secondary Structure Prediction: Chou–Fasman, Improved Chou–Fasman, and GOR

## 1. Project Description
This project implements three classical algorithms for protein secondary structure prediction:

- **Chou–Fasman**
A rule-based algorithm that uses predefined amino acid conformation propensities to identify helix, sheet, and coil segments via nucleation and extension rules.

- **Improved Chou–Fasman**
An updated version incorporating refined propensity values and more sophisticated segment detection thresholds to detect secondary structure.

- **GOR (Garnier–Osguthorpe–Robson)**
An information-theoretic method that evaluates the central residue and its surrounding residues to compute conditional probabilities for helix, sheet, and coil states.

All algorithms were implemented in **Go (Golang)**, and an **R Shiny** web application was developed to provide an interactive interface for sequence input, prediction execution, and visualization.

---

## 2. Folder Structure

project/
├── chou_fasman/
│   ├── Data/
│   │   ├── input_sequence.txt
│   ├── Tests/
│   ├── main.go
│   ├── datatype.go
│   ├── functions.go
│   ├── functions_test.go
│   ├── drawing.go
├── improved_chou_fasman/
│   ├── Tests/
│   ├── main.go
│   ├── data.go
│   ├── functions.go
│   ├── functions_test.go
│   ├── drawing.go
├── GOR/
│   ├── data/
│   │   ├── blindTest
│   │   ├── training
│   ├── test/
│   ├── train/
│   │   ├── main.go
│   ├── predict/
│   │   ├── drawing.go
│   │   ├── main.go
│   ├── Tests/
│   ├── gor.go
│   ├── functions_test.go
├── data/
├── app.R
└── README.md

> **Note:** The `Tests` folder for each algorithm contains testing data for subroutines.

---

## 3. Requirements

### Go Backend
- Go **version 1.24.5**
- Additional Go packages:
    - **Canvas** (package for drawing)
    Download from:
    https://programmingforlovers.com/wp-content/uploads/canvas.zip

### R Shiny Frontend
Install the following R packages:

```r
install.packages(c(
  "shiny",
  "shinyjs",
  "tidyverse"
))
```

## 4. Running GO code
### Chou-Fasman
```go
cd chou_fasman
go test // test for subroutines implemented in Chous-Fasman
./chou_fasman # run code
```

### Improved Chou–Fasman
```go
cd improved_chou_fasman
go test // test for subroutines implemented in CImproved Chou–Fasman
./improved_chou_fasman # run code
```

### GOR
```go
// Training the algorithm
cd GOR/train
go run train/main.go -ids data/training/list.txt -pssm_dir data/training/pssm -dssp_dir data/training/dssp -out gor_model.json


// Prediction (PSSM input)
go run ./predict \
  -model gor_model.json \
  -pssm "data/blindTest/pssm/6B8B:A.pssm"

// Prediction (Raw Sequence Input)
go run ./predict \
  -model gor_model.json \
  -seq "GSPRTVEEIFKDYSARRAALLRALTKDVDDFYSQCDPEKENLCLYGHPNESWEVNLPAEEVPPELPEPALGINFARDGMQRKDWLSLVAVHSDCWLLSVSFYFGARLNRNERKRLFSLINDLPTLFDVVTGRKAM"
```

## 5. RShiny
Run the Shiny app directly
Upload .txt file as input

Supported Input Types:
1. Single amino acid sequence (1 line)

Displays predicted structure for all three algorithms.

2. Sequence + experimental structure (2 line)

Shows predictions and accuracy comparison across the three algorithms.

3. Multiple entries: name + sequence + structure

Performs benchmarking comparison across multiple proteins.

## 6. Example for Results
### 1. Predicted secondary structure
![Predicted secondary structure](images/prediction_example.png)

### 2. Accuracy of the three prediction methods for one protein
![Accuracy comparison for one protein](images/accuracy_single_protein.png)

### 3. Accuracy of the three prediction methods across multiple proteins
![Accuracy comparison across multiple proteins](images/accuracy_multiple_proteins.png)


## 7. Notes
- Each algorithm has its own GO modules and tests
- Drawing functions require the Canvas package
- GOR requires running training step before prediction unless the trained parameters already exist.

## 8. Authors
This project was developed by:
- Bonnie Lai
- Kate Zhang
- Sorro Sun
- Yu-Lun Chen

## 9. Disclaimer
This README file has been edited for formatting and language clarity with the assistance of AI tools.
