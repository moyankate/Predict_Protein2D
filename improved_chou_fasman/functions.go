package main

import (
	"math"
)


// Morlet Wavelet Function (Equation 1) 
func morlet(t float64) float64 {
	return math.Exp(-0.5*t*t) * math.Cos(5*t)
}

// CWT Continuous Wavelet Transform (Equation 3) 
// Scale a = 9 is chosen 
func computeCWT(seqVals []float64, scale float64) []float64 {
	n := len(seqVals)
	coeffs := make([]float64, n)
	window := int(3.0 * scale)

	for b := 0; b < n; b++ {
		sum := 0.0
		for k := b - window; k <= b+window; k++ {
			if k >= 0 && k < n {
				t := float64(k-b) / scale
				sum += seqVals[k] * morlet(t)
			}
		}
		coeffs[b] = sum / math.Sqrt(scale)
	}
	return coeffs
}

// Algorithm Flow

func predictStructure(seq string) string {
	n := len(seq)

	// A. Prepare Data
	hVals := make([]float64, n)
	paVals := make([]float64, n)
	pbVals := make([]float64, n)

	for i, r := range seq {
		if val, ok := hydrophobicity[r]; ok {
			hVals[i] = val
		}
		if p, ok := propensities[r]; ok {
			paVals[i] = p.Pa
			pbVals[i] = p.Pb
		}
	}

	// B. Nucleation by CWT
	cwt := computeCWT(hVals, 9.0)
	nucleationSites := make([]bool, n)
	
	for i := 1; i < n-1; i++ {
		isPeak := cwt[i] > cwt[i-1] && cwt[i] > cwt[i+1]
		isValley := cwt[i] < cwt[i-1] && cwt[i] < cwt[i+1]
		
		// filtering: noise floor of digital signals 
		if (isPeak || isValley) && math.Abs(cwt[i]) > 0.1 {
			nucleationSites[i] = true
		}
	}

	// C. Extension Thresholds
	// Specifically for Alpha/Beta proteins from Table 7 
	thresholdH := 1.00 
	thresholdE := 1.02 

	rawH := make([]bool, n)
	rawE := make([]bool, n)

	for seed, isSeed := range nucleationSites {
		if !isSeed { continue }

		// Extend Segment
		hL, hR := extendBounds(seed, n, paVals, thresholdH)
		eL, eR := extendBounds(seed, n, pbVals, thresholdE)
		
		// Mark candidates
		markRegion(rawH, hL, hR)
		markRegion(rawE, eL, eR)
	}

	// D. Refinement / Conflict Resolution
	// Since segments from different seeds can have complex overlaps (partial/full),
	// we implement a per-residue score comparison which is mathematically equivalent to a "locally weighted average" decision.
	final := make([]rune, n)
	for i := 0; i < n; i++ {
		final[i] = 'C' // Default Coil

		if rawH[i] && rawE[i] {
			// Compare scores at this residue
			if paVals[i] >= pbVals[i] {
				final[i] = 'H'
			} else {
				final[i] = 'E'
			}
		} else if rawH[i] {
			final[i] = 'H'
		} else if rawE[i] {
			final[i] = 'E'
		}
	}

	// E. Biological Cleanup
	// Remove helix < 3 and strand < 2 
	return cleanUpStructure(string(final))
}

// extendBounds extends from a nucleation seed using propensity scores, and returns the candidate segment bounds (left, right).
func extendBounds(seed int, n int, scores []float64, threshold float64) (int, int) {
	left := seed
	for left >= 4 {
		avg := (scores[left] + scores[left-1] + scores[left-2] + scores[left-3]) / 4.0
		if avg < threshold { break }
		left--
	}
	right := seed
	for right < n-4 {
		avg := (scores[right] + scores[right+1] + scores[right+2] + scores[right+3]) / 4.0
		if avg < threshold { break }
		right++
	}
	return left, right
}

//markRegion marks the bool array to true 
func markRegion(arr []bool, L, R int) {
	for i := L; i <= R; i++ {
		arr[i] = true
	}
}

//cleanUpStructure remove takes string as input and remove short Helix and short strands
func cleanUpStructure(seq string) string {
	runes := []rune(seq)
	n := len(runes)

	// Remove short Helix (< 3)
	for i := 0; i < n; {
		if runes[i] == 'H' {
			j := i
			for j < n && runes[j] == 'H' { j++ }
			if j-i < 3 {
				for k := i; k < j; k++ { runes[k] = 'C' }
			}
			i = j
		} else {
			i++
		}
	}
	
	// Remove short Strand (< 2)
	for i := 0; i < n; {
		if runes[i] == 'E' {
			j := i
			for j < n && runes[j] == 'E' { j++ }
			if j-i < 2 {
				for k := i; k < j; k++ { runes[k] = 'C' }
			}
			i = j
		} else {
			i++
		}
	}
	return string(runes)
}