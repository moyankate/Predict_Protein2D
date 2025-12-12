// Author: Bonnie Lai
// Date: 2025-11-20
// Description: program to predict protein secondary structure using improved Chou-Fasman method

package main

import (
	"fmt"
	"os"
)

func main() {
	// test sequence
	if len(os.Args) < 2 {
		fmt.Println("Error: Please provide a file path")
		os.Exit(1)
	}
	filepath := os.Args[1]

	seq, err := ReadSequenceFromFile(filepath)

	if err != nil {
		fmt.Printf("Panic: %v\n", err)
		os.Exit(1)
	}

	fmt.Println("Improved Chou-Fasman Prediction (Paper-aligned)")
	fmt.Println("\nSequence:")
	fmt.Println(seq)
	fmt.Println("Sequence Length:", len(seq))

	result := PredictStructure(seq)

	fmt.Println("\nPrediction Result:")
	fmt.Println(result)

	fmt.Println("\nLegend:")
	fmt.Println("H = Alpha Helix")
	fmt.Println("E = Beta Strand (Sheet)")
	fmt.Println("C = Random Coil")

	// Draw bar chart
	err = drawPredictionBarChartCanvas(result, "improved_chou_fasman_prediction_bar_chart.png")
	if err != nil {
		fmt.Printf("Error drawing bar chart: %v\n", err)
		os.Exit(1)
	}
	fmt.Println("\nBar chart saved to improved_chou_fasman_prediction_bar_chart.png")
}
