// Author: Kate Zhang & Sorro Sun
// Date: 2025-11-20
// Description: program to predict protein secondary structure using GOR method
// Input: a PSSM file, a FASTA file, or a raw amino acid sequence string through command-line arguments
// Output: predicted secondary structure and a barchart image file

package main

import (
	"flag"
	"fmt"
	"log"
	"os"

	gor "protein_prediction/GOR"
)

func main() {
	modelFile := flag.String("model", "gor_model.json", "GOR model JSON file")
	pssmFile := flag.String("pssm", "", "PSSM file (optional, mutually exclusive with -fasta and -seq)")
	fastaFile := flag.String("fasta", "", "FASTA file (optional, mutually exclusive with -pssm and -seq)")
	seqString := flag.String("seq", "", "Directly type the amino acid sequence string (e.g. 'MVLSEGEWQL') (optional, mutually exclusive with -pssm and -fasta)")
	flag.Parse()

	if *pssmFile == "" && *fastaFile == "" && *seqString == "" {
		if len(os.Args) > 1 && os.Args[1] != "" {
			*fastaFile = os.Args[1]
		} else {
			log.Fatal("provide either -pssm OR -fasta OR -seq")
		}
	}

	model, err := gor.LoadGORModel(*modelFile)
	if err != nil {
		log.Fatalf("LoadGORModel error: %v", err)
	}

	var profile [][]float64

	if *pssmFile != "" {
		fmt.Printf("Loading PSSM from: %s\n", *pssmFile)
		profile, err = gor.ParsePSSM(*pssmFile)
		if err != nil {
			log.Fatalf("ParsePSSM error: %v", err)
		}
	} else if *fastaFile != "" {
		fmt.Printf("Loading FASTA from: %s\n", *fastaFile)
		seq, err := gor.ParseFASTA(*fastaFile)
		if err != nil {
			log.Fatalf("ParseFASTA error: %v", err)
		}
		profile = gor.SeqToProfile(seq)
	} else if *seqString != "" {
		fmt.Printf("Using raw sequence: %s\n", *seqString)
		profile = gor.SeqToProfile(*seqString)
	} else {
		log.Fatal("Please provide input via -pssm, -fasta, or -seq")
	}

	ss, err := gor.PredictGOR(model, profile)
	if err != nil {
		log.Fatalf("PredictGOR error: %v", err)
	}

	fmt.Println("Prediction Result:")
	fmt.Println(ss)

	// Canvas to draw barchart
	// PASS the predictedSequence INSTEAD of alphaFinal, betaFinal
	errC := drawPredictionBarChartCanvas(ss, "GOR_prediction_bar_chart.png")
	if errC != nil {
		fmt.Printf("Failed to draw a barchart: %v\n", errC)
	} else {
		fmt.Println("\n Successfully saved the barchart to 'GOR_prediction_bar_chart.png'")
	}

}
