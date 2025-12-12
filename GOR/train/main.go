package main

import (
	"flag"
	"fmt"
	"log"

	gor "GOR"
)

func main() {
	idListFile := flag.String("ids", "", "File with list of protein IDs")
	pssmDir := flag.String("pssm_dir", "", "Directory with <id>.pssm files")
	dsspDir := flag.String("dssp_dir", "", "Directory with <id>.dssp files")
	windowSize := flag.Int("window", 17, "Window size (default 17)")
	outModel := flag.String("out", "gor_model.json", "Output model file")
	flag.Parse()

	if *idListFile == "" || *pssmDir == "" || *dsspDir == "" {
		log.Fatal("please provide -ids, -pssm_dir, and -dssp_dir")
	}

	model, err := gor.TrainGOR(*idListFile, *pssmDir, *dsspDir, *windowSize)
	if err != nil {
		log.Fatalf("TrainGOR error: %v", err)
	}

	if err := gor.SaveGORModel(model, *outModel); err != nil {
		log.Fatalf("SaveGORModel error: %v", err)
	}

	fmt.Println("GOR model saved to", *outModel)
}
