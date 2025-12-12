package main

import (
	"flag"
	"fmt"
	"log"

	gor "GOR"
)

func main() {
	modelFile := flag.String("model", "gor_model.json", "GOR model JSON file")
	pssmFile := flag.String("pssm", "", "PSSM file (optional, mutually exclusive with -fasta)")
	fastaFile := flag.String("fasta", "", "FASTA file (optional, mutually exclusive with -pssm)")
	flag.Parse()

	if (*pssmFile == "" && *fastaFile == "") || (*pssmFile != "" && *fastaFile != "") {
		log.Fatal("provide either -pssm OR -fasta")
	}

	model, err := gor.LoadGORModel(*modelFile)
	if err != nil {
		log.Fatalf("LoadGORModel error: %v", err)
	}

	var profile [][]float64

	if *pssmFile != "" {
		profile, err = gor.ParsePSSM(*pssmFile)
		if err != nil {
			log.Fatalf("ParsePSSM error: %v", err)
		}
	} else {
		seq, err := gor.ParseFASTA(*fastaFile)
		if err != nil {
			log.Fatalf("ParseFASTA error: %v", err)
		}
		profile = gor.SeqToProfile(seq)
	}

	ss, err := gor.PredictGOR(model, profile)
	if err != nil {
		log.Fatalf("PredictGOR error: %v", err)
	}

	fmt.Println(ss)
}
