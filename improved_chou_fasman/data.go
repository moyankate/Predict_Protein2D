package main

// Table 1: Hydrophobicity values 
var hydrophobicity = map[rune]float64{
	'G': 0.00, 'Q': 0.00, 'S': 0.07, 'T': 0.07, 'N': 0.09,
	'D': 0.66, 'E': 0.67, 'R': 0.85, 'A': 0.87, 'H': 0.87,
	'C': 1.52, 'K': 1.64, 'M': 1.67, 'V': 1.87, 'L': 2.17,
	'Y': 2.76, 'P': 2.77, 'F': 2.87, 'I': 3.15, 'W': 3.77,
}

// Folding Type-Specific Propensities for Alpha/Beta Proteins
// From Jiang et al. (1998) Table I
type Propensity struct {
	Pa float64 // Alpha Helix propensity
	Pb float64 // Beta Strand propensity
}

var propensities = map[rune]Propensity{
	'A': {1.02, 0.83}, 'R': {0.98, 0.93}, 'N': {0.67, 0.89}, 'D': {1.01, 0.54},
	'C': {0.70, 1.19}, 'E': {1.51, 0.37}, 'Q': {1.11, 1.10}, 'G': {0.57, 0.75},
	'H': {1.00, 0.87}, 'I': {1.08, 1.60}, 'L': {1.21, 1.30}, 'K': {1.16, 0.74},
	'M': {1.45, 1.05}, 'F': {1.13, 1.38}, 'P': {0.57, 0.55}, 'S': {0.77, 0.75},
	'T': {0.83, 1.19}, 'W': {1.08, 1.37}, 'Y': {0.69, 1.47}, 'V': {1.06, 1.70},
}
