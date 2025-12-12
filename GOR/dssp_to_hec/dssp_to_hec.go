package main

import (
	"bufio"
	"fmt"
	"os"
	"strings"
)

// DSSPToHEC reads a DSSP file and returns a 3-state H/E/C sequence.
func DSSPToHEC(filename string) (string, error) {
	f, err := os.Open(filename)
	if err != nil {
		return "", err
	}
	defer f.Close()

	scanner := bufio.NewScanner(f)

	inTable := false
	var hec []rune

	for scanner.Scan() {
		line := scanner.Text()

		// Find the start of the residue table
		if !inTable {
			if strings.HasPrefix(line, "  #  RESIDUE AA STRUCTURE") {
				inTable = true
			}
			continue
		}

		// After the header line, all following lines are residue records
		// until EOF. Skip too-short lines just in case.
		if len(line) < 17 {
			continue
		}

		// DSSP structure code is in a fixed column (col 17 in classic format)
		// 0-based index 16.
		sym := rune(line[16])

		// Map DSSP 8-state code to 3-state H/E/C
		var c rune
		switch sym {
		// Helix-like
		case 'H', 'G', 'I':
			c = 'H'
		// Strand-like
		case 'E', 'B':
			c = 'E'
		// Everything else = coil
		default:
			c = 'C'
		}

		hec = append(hec, c)
	}

	if err := scanner.Err(); err != nil {
		return "", err
	}

	return string(hec), nil
}

func main() {
	if len(os.Args) < 2 {
		fmt.Println("Usage: dssp2hec <file.dssp>")
		os.Exit(1)
	}
	filename := os.Args[1]

	hec, err := DSSPToHEC(filename)
	if err != nil {
		fmt.Println("Error:", err)
		os.Exit(1)
	}

	fmt.Println(hec)
}
