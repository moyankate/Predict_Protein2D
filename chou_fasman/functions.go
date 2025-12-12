// Author: Yu-Lun Chen
// Date: 2025-12-11
// Description: Functions and subroutines for Chou-fasman.
// Comments are mainly created by AI tools

package main

import (
	"fmt"
	"math"
	"sort"
	"os"
	"strings"
)

// ConformationParameter retrieves the conformational parameter value for a specific amino acid 
// residue based on the target secondary structure.
func ConformationParameter(residue rune, conformation string) float64 {
	var paraMap map[rune]float64

	switch conformation {
	case "alpha":
		paraMap = alphaParam
	case "beta":
		paraMap = betaParam
	default:
		return 1.0
	}
	if val, exist := paraMap[residue]; exist {
		return val
	}
	return 1.0
}

// IsNucleationRegion checks if a given Region satisfies the criteria to be considered
// a potential nucleation site for a specific conformation (Alpha or Beta).
func IsNucleationRegion(region Region, conformation string, windowThreshold int, minParam float64) bool {
	seqRunes := []rune(SEQUENCE)
	window := seqRunes[region.Start:region.End]
	qualified := 0

	for _, residue := range window {
		if ConformationParameter(residue, conformation) > minParam {
			qualified++
		}
	}

	// The criteria include a minimum count of residues with high parameter values (windowThreshold)
	// and having the highest average parameter compared to other conformations.
	if qualified >= windowThreshold && ConformationHasHigestAverage(region, conformation) {
		return true
	}
	return false
}

// ConformationHasHigestAverage checks if the average conformational parameter
// for the given region and conformation is strictly greater than the average 
// parameters for all other major conformations.
func ConformationHasHigestAverage(region Region, conformation string) bool {
	avg := AverageParam(region, conformation)
	alphaAvg := AverageParam(region, "alpha")
	betaAvg := AverageParam(region, "beta")

	switch conformation {
	case "alpha":
		return avg >= betaAvg
	case "beta":
		return avg > alphaAvg
	default:
		return false
	}
}

// AverageParam calculates the average conformational parameter for a given 
// sequence region and conformation type.
func AverageParam(region Region, conformation string) float64 {
	seqRunes := []rune(SEQUENCE)
	if region.Start >= region.End || region.Start < 0 || region.End > len(seqRunes) {
		return 0.0
	}

	subSeq := seqRunes[region.Start:region.End]
	total := 0.0
	for _, residue := range subSeq {
		total += ConformationParameter(residue, conformation)
	}

	if len(subSeq) == 0 {
		return 0.0
	}
	return total / float64(len(subSeq))
}


// ----- Manage the nucleation region ----- //

// FindNucleationRegion iterates through the entire protein sequence using a sliding window
// of size 'windowSize' and identifies all sequence segments that qualify as a nucleation region
func FindNucleationRegion(windowThreshold, windowSize int, minParam float64, conformation string) []Region {
	var nucleationRegion []Region
	seqLen := len(SEQUENCE)
	i := 0
	j := windowSize

	for j <= seqLen {
		region := Region{Start: i, End: j}
		if IsNucleationRegion(region, conformation, windowThreshold, minParam) {
			nucleationRegion = append(nucleationRegion, region)
		}
		// shift window index
		i++
		j++
	}
	return nucleationRegion
}

// ExtendedRegions attempts to extend each found nucleation region bidirectionally 
// (N-terminal and C-terminal) as long as the 4-residue sliding window average 
// conformational parameter >= 1.0.
func ExtendedRegions(regions []Region, conformation string) []Region {
	var extendedRegions []Region
	for _, reg := range regions {
		N := 0 // delta index for N-terminal
		C := 0 // delta index for C-terminal

		// N-terminal extension: moves window towards Start=0
		shift := Region{Start: reg.Start, End: reg.Start + 4}
		for {
			newShift, ok := ShiftWindow(shift, 'N')
			if ! ok {
				break
			}
			if AverageParam(newShift, conformation) >= 1.0 {
				N++
				shift = newShift
			} else {
				break
			}
		}

		// C-terminal extension: moves window towards End=seqLen
		shift = Region{Start: reg.End - 4, End: reg.End}
		for {
			newShift, ok := ShiftWindow(shift, 'C')
			if ! ok {
				break
			}
			if AverageParam(newShift, conformation) >= 1.0 {
				C++
				shift = newShift
			} else {
				break
			}
		}
		extendedRegions = append(extendedRegions, Region{Start: reg.Start - N, End: reg.End + C})
	}
	return extendedRegions
}

// ShiftWindow moves a 4-residue window by one position towards the specified terminus 
// ('N' for N-terminal extension, 'C' for C-terminal extension) and checks boundary conditions.
func ShiftWindow(window Region, terminus rune) (Region, bool) {
	seqLen := len(SEQUENCE)

	if terminus == 'N' {
		if window.Start - 1 >= 0 {
			return Region{Start: window.Start - 1, End: window.End - 1}, true
		}
	}

	if terminus == 'C' {
		if window.End + 1 <= seqLen {
			return Region{Start: window.Start + 1, End: window.End + 1}, true
		}
	}

	// Returns the new window and a boolean indicating if the shift was successful.
	return Region{}, false
}


// ----- Filter extended region ---- //

// FilterExtendedRegions filters the extended regions, keeping only those whose
// overall average conformational parameter for the entire region >= the given threshold.
func FilterExtendedRegions(regions []Region, conformation string, threshold float64) []Region {
	if len(regions) == 0 {
		return []Region{}
	}

	var filteredRegions []Region
	for _, reg := range regions {
		if AverageParam(reg, conformation) >= threshold {
			filteredRegions = append(filteredRegions, reg)
		}
	}
	return filteredRegions
}

// ----- Merge functions ----- //

// mergeOverlappingRegions merges regions that overlap or are adjacent.
func MergeOverlappingRegions(regions []Region) []Region {
	if len(regions) < 2 {
		return regions
	}

	// sorted by region.Start for region in Regions
	sort.Slice(regions, func(i, j int) bool {
		return regions[i].Start < regions[j].Start
	})

	var mergedRegions []Region
	current := regions[0]

	for i := 1; i < len(regions); i++ {
		next := regions[i]
		if merged, ok := TryMerging(current, next); ok {
			current = merged
		} else {
			// combines overlapping regions into a single, larger region.
			mergedRegions = append(mergedRegions, current)
			current = next // begin with a new region
		}
	}
	mergedRegions = append(mergedRegions, current) // append the last region
	return mergedRegions
}

// tryMerging checks if two regions, i and j, overlap (i.End > j.Start).
func TryMerging(i, j Region) (Region, bool) {
	// if end_pos of i > start_pos of i+1: merge i & i+1: overlapped
	// returns a new merged region with the smallest Start and the largest End, along with 'true'.
	if i.End >= j.Start {
		// Start: the smallest start_pos, which is i.Start
		// End: max of i.end_pos & j.end_pos
		return Region{Start: int(math.Min(float64(i.Start), float64(j.Start))), End: int(math.Max(float64(i.End), float64(j.End)))}, true
	}
	// otherwise, it returns an empty Region and 'false'.
	return Region{}, false
}


// ----- Solving overlapped regions with different conformation ----- //

// solveOverlaps resolves conflicts between regions of two different conformations (conf1 vs conf2).
func SolveOverlaps(regions1 []Region, conf1 string, regions2 []Region, conf2 string) []Region {
	var solvedRegions []Region
	shouldKeep := make(map[Region]bool)

	for _, reg1 := range regions1 {
		shouldKeep[reg1] = true
	}

	// keeps a region from regions1 only if it is "stronger" (has higher average P) than any region from regions2 that it overlaps with.
	for _, reg1 := range regions1 {
		for _, reg2 := range regions2 {
			if AreOverlapping(reg1, conf1, reg2, conf2) {
				winner := StrongerAffinity(reg1, conf1, reg2, conf2)
				if winner == reg2 {
					// do not keep reg1
					shouldKeep[reg1] = false
					break
				}
			}
		}
	}

	for _, reg1 := range regions1 {
		if shouldKeep[reg1] == true {
			solvedRegions = append(solvedRegions, reg1)
		}
	}
	return solvedRegions
}

// AreOverlapping checks if two regions, regardless of their conformation, have any positional overlap.
func AreOverlapping(reg1 Region, conf1 string, reg2 Region, conf2 string) bool {
	return reg1.Start < reg2.End && reg2.Start < reg1.End
}

// strongerAffinity determines which of the two overlapping regions (reg1 or reg2)
// has a higher average conformational parameter across their entire length.
func StrongerAffinity(reg1 Region, conf1 string, reg2 Region, conf2 string) Region {
	avg1 := AverageParam(reg1, conf1)
	avg2 := AverageParam(reg2, conf2)

	if avg1 > avg2 {
		return reg1
	} else if avg2 > avg1 {
		return reg2
	} else {
		// if the averages are equal, Alpha Helix ('alpha') conformation wins.
		if conf1 == "alpha" {
			return reg1
		} else if conf2 == "alpha" {
			return reg2
		}
		// neither alpha nor beta for conf1/2
		return reg2
	}
}


// ----- Find Coils ----- //

// findCoilRegions identifies the contiguous regions that are neither Alpha nor Beta.
func FindCoilRegions(alphaFinal, betaFinal []Region) []Region {
    seqLen := len(SEQUENCE)
    if seqLen == 0 {
        return nil
    }
    
    // isPredicted tracks which indices are covered by Alpha or Beta
    isPredicted := make([]bool, seqLen)

    // Mark all Alpha residues as predicted
    for _, reg := range alphaFinal {
        for i := reg.Start; i < reg.End && i < seqLen; i++ {
            if i >= 0 {
                isPredicted[i] = true
            }
        }
    }

    // Mark all Beta residues as predicted
    for _, reg := range betaFinal {
        for i := reg.Start; i < reg.End && i < seqLen; i++ {
            if i >= 0 {
                isPredicted[i] = true
            }
        }
    }

    var coilRegions []Region
    inCoil := false
    currentStart := -1

    for i := 0; i < seqLen; i++ {
        if !isPredicted[i] {
            // Start of a new Coil region
            if !inCoil {
                inCoil = true
                currentStart = i
            }
        } else {
            // End of a Coil region
            if inCoil {
                coilRegions = append(coilRegions, Region{Start: currentStart, End: i})
                inCoil = false
                currentStart = -1
            }
        }
    }

    // Handle Coil region extending to the end of the sequence
    if inCoil {
        coilRegions = append(coilRegions, Region{Start: currentStart, End: seqLen})
    }

    return coilRegions
}


// ----- Produce sequence for secondary structure ----- //

// PredictSequenceFromRegions generates a secondary structure prediction sequence (H/E/C)
// based on the final resolved alpha helix and beta sheet regions.
// H = Helix (Alpha), E = Extended (Beta Sheet), C = Coil (Unpredicted)
func PredictSequenceFromRegions(alphaFinal, betaFinal []Region) string {
    seqLen := len(SEQUENCE)
    // Initialize prediction array with 'C' (Coil/Unpredicted)
    prediction := make([]rune, seqLen)
    for i := range prediction {
        prediction[i] = 'C'
    }

    // Step 1: Mark Beta regions as 'E'
    for _, reg := range betaFinal {
        for i := reg.Start; i < reg.End && i < seqLen; i++ {
            prediction[i] = 'E'
        }
    }
    
    // Step 2: Mark Alpha regions as 'H'
    for _, reg := range alphaFinal {
        for i := reg.Start; i < reg.End && i < seqLen; i++ {
            prediction[i] = 'H'
        }
    }

    return string(prediction)
}


// ----- Read data ----- //

// ReadSequenceFromFile reads the protein sequence from the specified file path.
func ReadSequenceFromFile(filepath string) (string, error) {
	content, err := os.ReadFile(filepath)
	if err != nil {
		return "", fmt.Errorf("Unable to read the file %s: %w", filepath, err)
	}
	
	rawSeq := string(content)
	cleanSeq := strings.Map(func(r rune) rune {
		if r == '\n' || r == '\r' || r == ' ' || r == '\t' {
			return -1
		}
		return r
	}, rawSeq)

	return strings.ToUpper(cleanSeq), nil
}
