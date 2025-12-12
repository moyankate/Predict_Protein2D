// Author: Yu-Lun Chen
// Date: 2025-12-11
// Description: Test functions for subroutines for Chou-fasman.
// The testing functions are partially derived from AI tools.
// Comments are mainly created by AI tools

package main

import (
	"bufio"
	"fmt"
	"math"
	"os"
	"sort"
	"strconv"
	"strings"
	"testing"
)

// NOTE: Global variable SEQUENCE is defined in datatype.go.
// We will temporarily overwrite it in tests for isolation and restore it afterwards.

// ----- General Helper Functions ----- //

// RegionsEqual checks if two slices of Region are equal by length and content.
func RegionsEqual(a, b []Region) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if a[i] != b[i] {
			return false
		}
	}
	return true
}

// ParseRegionString reads a comma-separated list of regions (e.g., "1:5, 10:12") into a slice of Region structs.
func ParseRegionString(s string) []Region {
	if s == "(empty)" || s == "" {
		return nil
	}
	var regions []Region
	pairs := strings.Split(s, ",")
	for _, pair := range pairs {
		parts := strings.Split(strings.TrimSpace(pair), ":")
		if len(parts) == 2 {
			start, err1 := strconv.Atoi(parts[0])
			end, err2 := strconv.Atoi(parts[1])
			if err1 == nil && err2 == nil {
				regions = append(regions, Region{Start: start, End: end})
			}
		}
	}
	return regions
}

// ReadTestFile loads the file content (lines starting with '# ') into a key-value map for structured parsing.
func ReadTestFile(filepath string) map[string]string {
	f, err := os.Open(filepath)
	if err != nil {
		panic(fmt.Sprintf("Error reading file %s: %v", filepath, err))
	}
	defer f.Close()

	data := make(map[string]string)
	scanner := bufio.NewScanner(f)
	
	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		if strings.HasPrefix(line, "#") {
			parts := strings.SplitN(strings.TrimPrefix(line, "# "), ":", 2)
			if len(parts) == 2 {
				key := strings.TrimSpace(parts[0])
				value := strings.TrimSpace(parts[1])
				data[key] = value
			}
		}
	}
	return data
}

// GetTestFiles reads all .txt files in the given directory and returns their full paths, sorted alphabetically.
func GetTestFiles(t *testing.T, dir string) []string {
	entries, err := os.ReadDir(dir)
	if err != nil {
		t.Fatalf("Failed to read test directory %s: %v. Please ensure the folder and files exist.", dir, err)
		return nil
	}
	var files []string
	for _, e := range entries {
		if !e.IsDir() && strings.HasSuffix(e.Name(), ".txt") {
			files = append(files, dir+"/"+e.Name())
		}
	}
	if len(files) == 0 {
		t.Logf("Warning: No .txt test files found in %s.", dir)
	}
	sort.Strings(files)
	return files
}


// --- 2. AverageParam ---

// AverageParamTest holds the structured input and expected output for the AverageParam function.
type AverageParamTest struct {
	ID           string
	Sequence     string
	Region       Region
	Conformation string
	Expected     float64
}

// ReadAverageParamTest reads a single test case for AverageParam from a file.
func ReadAverageParamTest(filepath string) AverageParamTest {
	data := ReadTestFile(filepath)

	tc := AverageParamTest{
		ID:           data["ID"],
		Sequence:     data["Sequence"],
		Conformation: data["Conformation"],
	}
	start, _ := strconv.Atoi(data["Start"])
	end, _ := strconv.Atoi(data["End"])
	tc.Region = Region{Start: start, End: end}
	tc.Expected, _ = strconv.ParseFloat(data["Expected"], 64)

	return tc
}

// TestAverageParam runs test cases for the AverageParam function from files in the specified directory.
func TestAverageParam(t *testing.T) {
	testDir := "Tests/averageParam"
	for _, filepath := range GetTestFiles(t, testDir) {
		tc := ReadAverageParamTest(filepath)
		t.Run(tc.ID, func(t *testing.T) {
			originalSequence := SEQUENCE
			SEQUENCE = tc.Sequence
			defer func() { SEQUENCE = originalSequence }()

			result := AverageParam(tc.Region, tc.Conformation)

			if math.Abs(result-tc.Expected) > 1e-6 {
				t.Errorf("Test %s failed.\nSequence: %s, Conf: %s.\nGot: %.4f, Expected: %.4f",
					tc.ID, tc.Sequence, tc.Conformation, result, tc.Expected)
			}
		})
	}
}

// --- 3. ConformationHasHigestAverage ---

// HighestAvgTest holds the structured input and expected output for the ConformationHasHigestAverage function.
type HighestAvgTest struct {
	ID           string
	Sequence     string
	Region       Region
	Conformation string
	Expected     bool
}

// ReadHighestAvgTest reads a single test case for ConformationHasHigestAverage from a file.
func ReadHighestAvgTest(filepath string) HighestAvgTest {
	data := ReadTestFile(filepath)
	
	tc := HighestAvgTest{
		ID:           data["ID"],
		Sequence:     data["Sequence"],
		Conformation: data["Conformation"],
		Expected:     data["Expected"] == "true",
	}
	start, _ := strconv.Atoi(data["Start"])
	end, _ := strconv.Atoi(data["End"])
	tc.Region = Region{Start: start, End: end}

	return tc
}

// TestConformationHasHigestAverage runs test cases for the ConformationHasHigestAverage function.
func TestConformationHasHigestAverage(t *testing.T) {
	testDir := "Tests/ConformationHasHigestAverage"
	for _, filepath := range GetTestFiles(t, testDir) {
		tc := ReadHighestAvgTest(filepath)
		t.Run(tc.ID, func(t *testing.T) {
			originalSequence := SEQUENCE
			SEQUENCE = tc.Sequence
			defer func() { SEQUENCE = originalSequence }()

			result := ConformationHasHigestAverage(tc.Region, tc.Conformation)

			if result != tc.Expected {
				t.Errorf("Test %s failed.\nSequence: %s, Conf: %s.\nGot: %v, Expected: %v",
					tc.ID, tc.Sequence, tc.Conformation, result, tc.Expected)
			}
		})
	}
}

// --- 4. IsNucleationRegion ---

// NucleationTest holds the structured input and expected output for the IsNucleationRegion function.
type NucleationTest struct {
	ID              string
	Sequence        string
	Region          Region
	Conformation    string
	WindowThreshold int
	MinParam        float64
	Expected        bool
}

// ReadNucleationTest reads a single test case for IsNucleationRegion from a file.
func ReadNucleationTest(filepath string) NucleationTest {
	data := ReadTestFile(filepath)

	tc := NucleationTest{
		ID:           data["ID"],
		Sequence:     data["Sequence"],
		Conformation: data["Conformation"],
		Expected:     data["Expected"] == "true",
	}
	start, _ := strconv.Atoi(data["Start"])
	end, _ := strconv.Atoi(data["End"])
	tc.Region = Region{Start: start, End: end}
	tc.WindowThreshold, _ = strconv.Atoi(data["WindowThreshold"])
	tc.MinParam, _ = strconv.ParseFloat(data["MinParam"], 64)

	return tc
}

// TestIsNucleationRegion runs test cases for the IsNucleationRegion function.
func TestIsNucleationRegion(t *testing.T) {
	testDir := "Tests/IsNucleationRegion"
	for _, filepath := range GetTestFiles(t, testDir) {
		tc := ReadNucleationTest(filepath)
		t.Run(tc.ID, func(t *testing.T) {
			originalSequence := SEQUENCE
			SEQUENCE = tc.Sequence
			defer func() { SEQUENCE = originalSequence }()

			result := IsNucleationRegion(tc.Region, tc.Conformation, tc.WindowThreshold, tc.MinParam)

			if result != tc.Expected {
				t.Errorf("Test %s failed.\nGot: %v, Expected: %v",
							tc.ID, result, tc.Expected)
			}
		})
	}
}

// --- 5. FindNucleationRegion ---

// FindNucleationTest holds the structured input and expected output for the FindNucleationRegion function.
type FindNucleationTest struct {
	ID              string
	Sequence        string
	Conformation    string
	WindowThreshold int
	WindowSize      int
	MinParam        float64
	Expected        []Region
}

// ReadFindNucleationTest reads a single test case for FindNucleationRegion from a file.
func ReadFindNucleationTest(filepath string) FindNucleationTest {
	data := ReadTestFile(filepath)

	tc := FindNucleationTest{
		ID:           data["ID"],
		Sequence:     data["Sequence"],
		Conformation: data["Conformation"],
	}
	tc.WindowThreshold, _ = strconv.Atoi(data["WindowThreshold"])
	tc.WindowSize, _ = strconv.Atoi(data["WindowSize"])
	tc.MinParam, _ = strconv.ParseFloat(data["MinParam"], 64)
	tc.Expected = ParseRegionString(data["ExpectedRegions"])

	return tc
}

// TestFindNucleationRegion runs test cases for the FindNucleationRegion function.
func TestFindNucleationRegion(t *testing.T) {
	testDir := "Tests/FindNucleationRegion"
	for _, filepath := range GetTestFiles(t, testDir) {
		tc := ReadFindNucleationTest(filepath)
		t.Run(tc.ID, func(t *testing.T) {
			originalSequence := SEQUENCE
			SEQUENCE = tc.Sequence
			defer func() { SEQUENCE = originalSequence }()

			result := FindNucleationRegion(tc.WindowThreshold, tc.WindowSize, tc.MinParam, tc.Conformation)

			if !RegionsEqual(result, tc.Expected) {
				t.Errorf("Test %s failed.\nGot: %v, Expected: %v",
					tc.ID, result, tc.Expected)
			}
		})
	}
}

// --- 6. ShiftWindow ---

// ShiftWindowTest holds the structured input and expected output for the ShiftWindow function.
type ShiftWindowTest struct {
	ID           string
	Sequence     string // Needed for sequence length boundary check
	Window       Region
	Terminus     rune
	Expected     Region
	ExpectedOK   bool
}

// ReadShiftWindowTest reads a single test case for ShiftWindow from a file.
func ReadShiftWindowTest(filepath string) ShiftWindowTest {
	data := ReadTestFile(filepath)

	tc := ShiftWindowTest{
		ID:       data["ID"],
		Sequence: data["Sequence"],
		Terminus: []rune(data["Terminus"])[0],
		ExpectedOK: data["ExpectedOK"] == "true",
	}
	start, _ := strconv.Atoi(data["Start"])
	end, _ := strconv.Atoi(data["End"])
	tc.Window = Region{Start: start, End: end}

	expectedStart, _ := strconv.Atoi(data["ExpectedStart"])
	expectedEnd, _ := strconv.Atoi(data["ExpectedEnd"])
	if tc.ExpectedOK {
		tc.Expected = Region{Start: expectedStart, End: expectedEnd}
	} else {
		tc.Expected = Region{}
	}

	return tc
}

// TestShiftWindow runs test cases for the ShiftWindow function.
func TestShiftWindow(t *testing.T) {
	testDir := "Tests/ShiftWindow"
	for _, filepath := range GetTestFiles(t, testDir) {
		tc := ReadShiftWindowTest(filepath)
		t.Run(tc.ID, func(t *testing.T) {
			originalSequence := SEQUENCE
			SEQUENCE = tc.Sequence
			defer func() { SEQUENCE = originalSequence }()

			result, ok := ShiftWindow(tc.Window, tc.Terminus)

			if result != tc.Expected || ok != tc.ExpectedOK {
				t.Errorf("Test %s failed. Terminus: %c\nGot: %v, %v\nExpected: %v, %v",
					tc.ID, tc.Terminus, result, ok, tc.Expected, tc.ExpectedOK)
			}
		})
	}
}

// --- 7. ExtendedRegions ---

// ExtendedRegionsTest holds the structured input and expected output for the ExtendedRegions function.
type ExtendedRegionsTest struct {
	ID           string
	Sequence     string
	Regions      []Region
	Conformation string
	Expected     []Region
}

// ReadExtendedRegionsTest reads a single test case for ExtendedRegions from a file.
func ReadExtendedRegionsTest(filepath string) ExtendedRegionsTest {
	data := ReadTestFile(filepath)

	tc := ExtendedRegionsTest{
		ID:           data["ID"],
		Sequence:     data["Sequence"],
		Conformation: data["Conformation"],
	}
	tc.Regions = ParseRegionString(data["InputRegions"])
	tc.Expected = ParseRegionString(data["ExpectedRegions"])

	return tc
}

// TestExtendedRegions runs test cases for the ExtendedRegions function.
func TestExtendedRegions(t *testing.T) {
	testDir := "Tests/ExtendedRegions"
	for _, filepath := range GetTestFiles(t, testDir) {
		tc := ReadExtendedRegionsTest(filepath)
		t.Run(tc.ID, func(t *testing.T) {
			originalSequence := SEQUENCE
			SEQUENCE = tc.Sequence
			defer func() { SEQUENCE = originalSequence }()

			result := ExtendedRegions(tc.Regions, tc.Conformation)

			if !RegionsEqual(result, tc.Expected) {
				t.Errorf("Test %s failed. Got: %v, Expected: %v",
					tc.ID, result, tc.Expected)
			}
		})
	}
}


// --- 8. FilterExtendedRegions ---

// FilterExtendedTest holds the structured input and expected output for the FilterExtendedRegions function.
type FilterExtendedTest struct {
	ID           string
	Sequence     string // Needed for averageParam calls within the function
	Regions      []Region
	Conformation string
	Threshold    float64
	Expected     []Region
}

// ReadFilterExtendedTest reads a single test case for FilterExtendedRegions from a file.
func ReadFilterExtendedTest(filepath string) FilterExtendedTest {
	data := ReadTestFile(filepath)
	
	tc := FilterExtendedTest{
		ID:           data["ID"],
		Sequence:     data["Sequence"],
		Conformation: data["Conformation"],
	}
	tc.Regions = ParseRegionString(data["InputRegions"])
	tc.Threshold, _ = strconv.ParseFloat(data["Threshold"], 64)
	tc.Expected = ParseRegionString(data["ExpectedRegions"])

	return tc
}

// TestFilterExtendedRegions runs test cases for the FilterExtendedRegions function.
func TestFilterExtendedRegions(t *testing.T) {
	testDir := "Tests/FilterExtendedRegions"
	for _, filepath := range GetTestFiles(t, testDir) {
		tc := ReadFilterExtendedTest(filepath)
		t.Run(tc.ID, func(t *testing.T) {
			originalSequence := SEQUENCE
			SEQUENCE = tc.Sequence
			defer func() { SEQUENCE = originalSequence }()

			result := FilterExtendedRegions(tc.Regions, tc.Conformation, tc.Threshold)

			if !RegionsEqual(result, tc.Expected) {
				t.Errorf("Test %s failed. Got: %v, Expected: %v",
					tc.ID, result, tc.Expected)
			}
		})
	}
}

// --- 9. TryMerging ---

// TryMergingTest holds the structured input and expected output for the TryMerging function.
type TryMergingTest struct {
	ID           string
	Region1      Region
	Region2      Region
	Expected     Region
	ExpectedOK   bool
}

// ReadTryMergingTest reads a single test case for TryMerging from a file.
func ReadTryMergingTest(filepath string) TryMergingTest {
	data := ReadTestFile(filepath)

	tc := TryMergingTest{
		ID:       data["ID"],
		ExpectedOK: data["ExpectedOK"] == "true",
	}
	start1, _ := strconv.Atoi(data["Start1"])
	end1, _ := strconv.Atoi(data["End1"])
	tc.Region1 = Region{Start: start1, End: end1}
	start2, _ := strconv.Atoi(data["Start2"])
	end2, _ := strconv.Atoi(data["End2"])
	tc.Region2 = Region{Start: start2, End: end2}

	expectedStart, _ := strconv.Atoi(data["ExpectedStart"])
	expectedEnd, _ := strconv.Atoi(data["ExpectedEnd"])
	if tc.ExpectedOK {
		tc.Expected = Region{Start: expectedStart, End: expectedEnd}
	} else {
		tc.Expected = Region{}
	}

	return tc
}

// TestTryMerging runs test cases for the TryMerging function.
func TestTryMerging(t *testing.T) {
	testDir := "Tests/TryMerging"
	for _, filepath := range GetTestFiles(t, testDir) {
		tc := ReadTryMergingTest(filepath)
		t.Run(tc.ID, func(t *testing.T) {
			result, ok := TryMerging(tc.Region1, tc.Region2)

			if result != tc.Expected || ok != tc.ExpectedOK {
				t.Errorf("Test %s failed. R1: %v, R2: %v\nGot: %v, %v\nExpected: %v, %v",
					tc.ID, tc.Region1, tc.Region2, result, ok, tc.Expected, tc.ExpectedOK)
			}
		})
	}
}


// --- 10. MergeOverlappingRegions ---

// MergeOverlapTest holds the structured input and expected output for the MergeOverlappingRegions function.
type MergeOverlapTest struct {
	ID             string
	InputRegions   []Region
	Expected       []Region
}

// ReadMergeOverlapTest reads a single test case for MergeOverlappingRegions from a file.
func ReadMergeOverlapTest(filepath string) MergeOverlapTest {
	data := ReadTestFile(filepath)

	tc := MergeOverlapTest{
		ID: data["ID"],
	}
	tc.InputRegions = ParseRegionString(data["InputRegions"])
	tc.Expected = ParseRegionString(data["ExpectedRegions"])

	return tc
}

// TestMergeOverlappingRegions runs test cases for the MergeOverlappingRegions function.
func TestMergeOverlappingRegions(t *testing.T) {
	testDir := "Tests/MergeOverlappingRegions"
	for _, filepath := range GetTestFiles(t, testDir) {
		tc := ReadMergeOverlapTest(filepath)
		t.Run(tc.ID, func(t *testing.T) {
			result := MergeOverlappingRegions(tc.InputRegions)

			if !RegionsEqual(result, tc.Expected) {
				t.Errorf("Test %s failed. Got: %v\nExpected: %v",
					tc.ID, result, tc.Expected)
			}
		})
	}
}


// --- 11. AreOverlapping ---

// AreOverlappingTest holds the structured input and expected output for the AreOverlapping function.
type AreOverlappingTest struct {
	ID      string
	Region1 Region
	Conf1   string // Not used in function logic, but required by signature
	Region2 Region
	Conf2   string // Not used in function logic, but required by signature
	Expected bool
}

// ReadAreOverlappingTest reads a single test case for AreOverlapping from a file.
func ReadAreOverlappingTest(filepath string) AreOverlappingTest {
	data := ReadTestFile(filepath)

	tc := AreOverlappingTest{
		ID:       data["ID"],
		Conf1:    data["Conf1"],
		Conf2:    data["Conf2"],
		Expected: data["Expected"] == "true",
	}
	start1, _ := strconv.Atoi(data["Region1Start"])
	end1, _ := strconv.Atoi(data["Region1End"])
	tc.Region1 = Region{Start: start1, End: end1}
	start2, _ := strconv.Atoi(data["Region2Start"])
	end2, _ := strconv.Atoi(data["Region2End"])
	tc.Region2 = Region{Start: start2, End: end2}

	return tc
}

// TestAreOverlapping runs test cases for the AreOverlapping function.
func TestAreOverlapping(t *testing.T) {
	testDir := "Tests/AreOverlapping"
	for _, filepath := range GetTestFiles(t, testDir) {
		tc := ReadAreOverlappingTest(filepath)
		t.Run(tc.ID, func(t *testing.T) {
			result := AreOverlapping(tc.Region1, tc.Conf1, tc.Region2, tc.Conf2)

			if result != tc.Expected {
				t.Errorf("Test %s failed. R1: %v, R2: %v\nGot: %v, Expected: %v",
					tc.ID, tc.Region1, tc.Region2, result, tc.Expected)
			}
		})
	}
}

// --- 12. StrongerAffinity ---

// StrongerAffinityTest holds the structured input and expected output for the StrongerAffinity function.
type StrongerAffinityTest struct {
	ID      string
	Sequence string
	Region1 Region
	Conf1   string
	Region2 Region
	Conf2   string
	Expected Region // Expected winner region (must match R1 or R2)
}

// ReadStrongerAffinityTest reads a single test case for StrongerAffinity from a file.
func ReadStrongerAffinityTest(filepath string) StrongerAffinityTest {
    data := ReadTestFile(filepath)

    tc := StrongerAffinityTest{
        ID:       data["ID"],
        Sequence: data["Sequence"],
        Conf1:    data["Conf1"],
        Conf2:    data["Conf2"],
    }

    start1, _ := strconv.Atoi(data["Region1Start"])
    end1, _ := strconv.Atoi(data["Region1End"])
    tc.Region1 = Region{Start: start1, End: end1}

    start2, _ := strconv.Atoi(data["Region2Start"])
    end2, _ := strconv.Atoi(data["Region2End"])
    tc.Region2 = Region{Start: start2, End: end2}

    expectedID := data["ExpectedWinner"]
	if expectedID == "R1" {
		tc.Expected = tc.Region1
	} else if expectedID == "R2" {
		tc.Expected = tc.Region2
	} else {
		// Should not happen for a well-formed test, but handle potential error
		panic(fmt.Errorf("Invalid ExpectedWinner tag '%s' in test file %s. Must be 'R1' or 'R2'.", expectedID, filepath))
	}

	return tc
}

// TestStrongerAffinity runs test cases for the StrongerAffinity function.
func TestStrongerAffinity(t *testing.T) {
	testDir := "Tests/StrongerAffinity"
	for _, filepath := range GetTestFiles(t, testDir) {
		tc := ReadStrongerAffinityTest(filepath)
		t.Run(tc.ID, func(t *testing.T) {
			originalSequence := SEQUENCE
			SEQUENCE = tc.Sequence
			defer func() { SEQUENCE = originalSequence }()

			result := StrongerAffinity(tc.Region1, tc.Conf1, tc.Region2, tc.Conf2)

			if result != tc.Expected {
				t.Errorf("Test %s failed.\nConf1: %s, Conf2: %s\nGot: %v, Expected: %v",
					tc.ID, tc.Conf1, tc.Conf2, result, tc.Expected)
			}
		})
	}
}

// --- 13. solveOverlaps ---

// SolveOverlapsTest holds the structured input and expected output for the SolveOverlaps function.
type SolveOverlapsTest struct {
	ID       string
	Sequence string
	Regions1 []Region // Conf1 regions to keep/filter
	Conf1    string
	Regions2 []Region // Conf2 regions (competitors)
	Conf2    string
	Expected []Region // Expected remaining regions in Regions1
}

// ReadSolveOverlapsTest reads a single test case for SolveOverlaps from a file.
func ReadSolveOverlapsTest(filepath string) SolveOverlapsTest {
	data := ReadTestFile(filepath)

	tc := SolveOverlapsTest{
		ID:       data["ID"],
		Sequence: data["Sequence"],
		Conf1:    data["Conf1"],
		Conf2:    data["Conf2"],
	}
	tc.Regions1 = ParseRegionString(data["Regions1"])
	tc.Regions2 = ParseRegionString(data["Regions2"])
	tc.Expected = ParseRegionString(data["ExpectedRegions"])

	return tc
}

// TestSolveOverlaps runs test cases for the SolveOverlaps function.
func TestSolveOverlaps(t *testing.T) {
	testDir := "Tests/SolveOverlaps"
	for _, filepath := range GetTestFiles(t, testDir) {
		tc := ReadSolveOverlapsTest(filepath)
		t.Run(tc.ID, func(t *testing.T) {
			originalSequence := SEQUENCE
			SEQUENCE = tc.Sequence
			defer func() { SEQUENCE = originalSequence }()

			result := SolveOverlaps(tc.Regions1, tc.Conf1, tc.Regions2, tc.Conf2)

			if !RegionsEqual(result, tc.Expected) {
				t.Errorf("Test %s failed.\nConf1: %s, Conf2: %s\nGot: %v\nExpected: %v",
					tc.ID, tc.Conf1, tc.Conf2, result, tc.Expected)
			}
		})
	}
}

// --- 14. FindCoilRegions ---

// FindCoilTest holds the structured input and expected output for the FindCoilRegions function.
type FindCoilTest struct {
	ID          string
	Sequence    string
	AlphaFinal  []Region
	BetaFinal   []Region
	Expected    []Region
}

// ReadFindCoilTest reads a single test case for FindCoilRegions from a file.
func ReadFindCoilTest(filepath string) FindCoilTest {
	data := ReadTestFile(filepath)

	tc := FindCoilTest{
		ID:       data["ID"],
		Sequence: data["Sequence"],
	}
	tc.AlphaFinal = ParseRegionString(data["AlphaFinal"])
	tc.BetaFinal = ParseRegionString(data["BetaFinal"])
	tc.Expected = ParseRegionString(data["ExpectedRegions"])

	return tc
}

// TestFindCoilRegions runs test cases for the FindCoilRegions function.
func TestFindCoilRegions(t *testing.T) {
	testDir := "Tests/FindCoilRegions"
	for _, filepath := range GetTestFiles(t, testDir) {
		tc := ReadFindCoilTest(filepath)
		t.Run(tc.ID, func(t *testing.T) {
			originalSequence := SEQUENCE
			SEQUENCE = tc.Sequence
			defer func() { SEQUENCE = originalSequence }()

			result := FindCoilRegions(tc.AlphaFinal, tc.BetaFinal)

			if !RegionsEqual(result, tc.Expected) {
				t.Errorf("Test %s failed. Sequence Length: %d\nGot: %v\nExpected: %v",
					tc.ID, len(tc.Sequence), result, tc.Expected)
			}
		})
	}
}

// --- 15. predictSequenceFromRegions ---

// PredictSequenceTest holds the structured input and expected output for the PredictSequenceFromRegions function.
type PredictSequenceTest struct {
	ID          string
	Sequence    string
	AlphaFinal  []Region
	BetaFinal   []Region
	Expected    string
}

// ReadPredictSequenceTest reads a single test case for PredictSequenceFromRegions from a file.
func ReadPredictSequenceTest(filepath string) PredictSequenceTest {
	data := ReadTestFile(filepath)

	tc := PredictSequenceTest{
		ID:       data["ID"],
		Sequence: data["Sequence"],
		Expected: data["ExpectedPrediction"],
	}
	tc.AlphaFinal = ParseRegionString(data["AlphaFinal"])
	tc.BetaFinal = ParseRegionString(data["BetaFinal"])

	return tc
}

// TestPredictSequenceFromRegions runs test cases for the PredictSequenceFromRegions function.
func TestPredictSequenceFromRegions(t *testing.T) {
	testDir := "Tests/PredictSequenceFromRegions"
	for _, filepath := range GetTestFiles(t, testDir) {
		tc := ReadPredictSequenceTest(filepath)
		t.Run(tc.ID, func(t *testing.T) {
			originalSequence := SEQUENCE
			SEQUENCE = tc.Sequence
			defer func() { SEQUENCE = originalSequence }()

			result := PredictSequenceFromRegions(tc.AlphaFinal, tc.BetaFinal)

			if result != tc.Expected {
				t.Errorf("Test %s failed.\nSequence: %s\nGot: %s\nExpected: %s",
					tc.ID, tc.Sequence, result, tc.Expected)
			}
		})
	}
}
