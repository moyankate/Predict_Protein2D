package main

import (
	"fmt"
)


func main() {
	seq := "MVLSEGEWQLVLHVWAKVEADVAGHGQDILIRLFKSHPETLEKFDRVKHLKTEAEMKASEDLKKHGVTVLTALGAILKKKGHHEAELKPLAQSHATKHKIPIKYLEFISEAIIHVLHSRHPGDFGADAQGAMNKALELFRKDMASNYKELGFQG"

	fmt.Println("Improved Chou-Fasman Prediction (Paper-aligned)")
	fmt.Println("\nSequence:")
	fmt.Println(seq)
	fmt.Println("Sequence Length:", len(seq))

	result := predictStructure(seq)

	fmt.Println("\nPrediction Result:")
	fmt.Println(result)

	fmt.Println("\nLegend:")
	fmt.Println("H = Alpha Helix")
	fmt.Println("E = Beta Strand (Sheet)")
	fmt.Println("C = Random Coil")
}



