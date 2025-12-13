// Author: Yu-Lun Chen
// Date: 2025-12-11
// Description: Functions and subroutines for Chou-fasman.
// Comments are mainly created by AI tools

package main

import (
	"fmt"
	"os"
)

// parameter
var alphaParam = map[rune]float64{
	'A': 1.42, 'R': 0.98, 'N': 0.67, 'D': 1.01, 'C': 0.70,
	'Q': 1.11, 'E': 1.51, 'G': 0.57, 'H': 1.00, 'I': 1.08,
	'L': 1.21, 'K': 1.16, 'M': 1.45, 'F': 1.13, 'P': 0.57,
	'S': 0.77, 'T': 0.83, 'W': 1.08, 'Y': 0.69, 'V': 1.06,
	'X': 1.00,
}

var betaParam = map[rune]float64{
	'A': 0.83, 'R': 0.93, 'N': 0.89, 'D': 0.54, 'C': 1.19,
	'Q': 1.10, 'E': 0.37, 'G': 0.75, 'H': 0.87, 'I': 1.60,
	'L': 1.30, 'K': 0.74, 'M': 1.05, 'F': 1.38, 'P': 0.55,
	'S': 0.75, 'T': 1.19, 'W': 1.37, 'Y': 1.47, 'V': 1.70,
	'X': 1.00,
}

func main() {
	// test sequence is stored in "Data/fileName.txt"
	const filepath = "Data/5Y53B.txt"
	seq, err := ReadSequenceFromFile(filepath)
	if err != nil {
		fmt.Printf("Error: %v\n", err)
		os.Exit(1)
	}
	SEQUENCE = seq

	// 1. calculate the possible nucleation region
	// Alpha: window size = 6, threshold = 4, min param = 1
	alphaRegions := FindNucleationRegion(4, 6, 1.0, "alpha")
	// Beta: window size = 5, threshold = 3, min param = 1
	betaRegions := FindNucleationRegion(3, 5, 1.0, "beta")

	// 2. Extend each nucleation
	alphaRegions = ExtendedRegions(alphaRegions, "alpha")
	betaRegions = ExtendedRegions(betaRegions, "beta")

	// 3. Filter out the extended regions
	alphaRegions = FilterExtendedRegions(alphaRegions, "alpha", 1.03)
	betaRegions = FilterExtendedRegions(betaRegions, "beta", 1.05)

	// 4. Merge extended regions
	alphaRegions = MergeOverlappingRegions(alphaRegions)
	betaRegions = MergeOverlappingRegions(betaRegions)

	// 5. Solve the overlapped regions
	alphaFinal := SolveOverlaps(alphaRegions, "alpha", betaRegions, "beta")
	betaFinal := SolveOverlaps(betaRegions, "beta", alphaRegions, "alpha")

	// 6. Find Coil regions (new step)
	coilFinal := FindCoilRegions(alphaFinal, betaFinal)

	// 7. Generate the full prediction sequence (H/E/C)
	predictedSequence := PredictSequenceFromRegions(alphaFinal, betaFinal)

	// ----- Output ----- //
	fmt.Println("--- Final Output ---")

	fmt.Println("Alpha regions:")
	fmt.Println("Alpha (H): ", alphaFinal)

	fmt.Println("\nBeta regions:")
	fmt.Println("Beta (E): ", betaFinal)

	fmt.Println("\nCoil regions:")
	fmt.Println("Coil (C): ", coilFinal)

	fmt.Println("\nPredicted Sequence (H/E/C):")
	fmt.Println(predictedSequence)

	// Canvas to draw barchart
	// PASS the predictedSequence INSTEAD of alphaFinal, betaFinal
	errC := drawPredictionBarChartCanvas(predictedSequence, "chou_fasman_prediction_bar_chart.png")
	if errC != nil {
		fmt.Printf("Failed to draw a barchart: %v\n", errC)
	} else {
		fmt.Println("\nSuccessfully saved the barchart to 'chou_fasman_prediction_bar_chart.png'")
	}
}
