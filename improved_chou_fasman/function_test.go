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

// --- General Helper Functions ---
// round6 rounds a float to 6 decimal places.
func round6(x float64) float64 {
	return math.Round(x*1e6) / 1e6
}

// ParseFloatSlice reads a comma-separated list of floats (e.g., "0.5, 1.2, -3.4")
func ParseFloatSlice(s string) []float64 {
	if s == "(empty)" || s == "" {
		return nil
	}
	var floats []float64
	parts := strings.Split(s, ",")
	for _, p := range parts {
		f, err := strconv.ParseFloat(strings.TrimSpace(p), 64)
		if err == nil {
			floats = append(floats, f)
		}
	}
	return floats
}

// ReadTestFile loads the file content into a map for structured parsing.
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

// GetTestFiles reads all .txt files in the given directory.
func GetTestFiles(t *testing.T, dir string) []string {
	// Check if directory exists, if not create it just to avoid panic during initial run if empty
	if _, err := os.Stat(dir); os.IsNotExist(err) {
		t.Logf("Directory %s does not exist, skipping tests in this dir.", dir)
		return nil
	}

	entries, err := os.ReadDir(dir)
	if err != nil {
		t.Fatalf("Failed to read test directory %s: %v", dir, err)
		return nil
	}
	var files []string
	for _, e := range entries {
		if !e.IsDir() && strings.HasSuffix(e.Name(), ".txt") {
			files = append(files, dir+"/"+e.Name())
		}
	}
	sort.Strings(files)
	return files
}

// --- 1. Morlet ---

type MorletTest struct {
	ID       string
	T        float64
	Expected float64
}

func ReadMorletTest(filepath string) MorletTest {
	data := ReadTestFile(filepath)
	tc := MorletTest{
		ID: data["ID"],
	}
	tc.T, _ = strconv.ParseFloat(data["T"], 64)
	tc.Expected, _ = strconv.ParseFloat(data["Expected"], 64)
	return tc
}

func TestMorlet(t *testing.T) {
	testDir := "Tests/Morlet"
	for _, filepath := range GetTestFiles(t, testDir) {
		tc := ReadMorletTest(filepath)
		t.Run(tc.ID, func(t *testing.T) {
			result := Morlet(tc.T)
			
			if round6(result) != round6(tc.Expected) {
				t.Errorf("Test %s failed.\nInput T: %.4f\nGot: %.6f, Expected: %.6f",
					tc.ID, tc.T, result, tc.Expected)
			}
		})
	}
}

// --- 2. ComputeCWT ---

type ComputeCWTTest struct {
	ID       string
	SeqVals  []float64
	Scale    float64
	Expected []float64
}

func ReadComputeCWTTest(filepath string) ComputeCWTTest {
	data := ReadTestFile(filepath)
	tc := ComputeCWTTest{
		ID: data["ID"],
	}
	tc.SeqVals = ParseFloatSlice(data["InputVals"])
	tc.Scale, _ = strconv.ParseFloat(data["Scale"], 64)
	tc.Expected = ParseFloatSlice(data["ExpectedVals"])
	return tc
}

func TestComputeCWT(t *testing.T) {
	testDir := "Tests/ComputeCWT"
	for _, filepath := range GetTestFiles(t, testDir) {
		tc := ReadComputeCWTTest(filepath)
		t.Run(tc.ID, func(t *testing.T) {
			result := ComputeCWT(tc.SeqVals, tc.Scale)

			// Debug print if lengths differ
			if len(result) != len(tc.Expected) {
				t.Fatalf("Test %s failed length check. Got len %d, expected %d", tc.ID, len(result), len(tc.Expected))
			}

			// Check the value
			for i := range result {
				if round6(result[i]) != round6(tc.Expected[i]) {
					t.Errorf("Test %s failed at index %d.\nGot:%.6f\nExpected: %.6f", 
						tc.ID, i, result[i], tc.Expected[i])
				}
			}
		})
	}
}

// --- 3. ExtendBounds ---

type ExtendBoundsTest struct {
	ID        string
	Seed      int
	N         int
	Scores    []float64
	Threshold float64
	ExpectedL int
	ExpectedR int
}

func ReadExtendBoundsTest(filepath string) ExtendBoundsTest {
	data := ReadTestFile(filepath)
	tc := ExtendBoundsTest{
		ID: data["ID"],
	}
	tc.Seed, _ = strconv.Atoi(data["Seed"])
	tc.N, _ = strconv.Atoi(data["N"])
	tc.Scores = ParseFloatSlice(data["Scores"])
	tc.Threshold, _ = strconv.ParseFloat(data["Threshold"], 64)
	tc.ExpectedL, _ = strconv.Atoi(data["ExpectedL"])
	tc.ExpectedR, _ = strconv.Atoi(data["ExpectedR"])
	return tc
}

func TestExtendBounds(t *testing.T) {
	testDir := "Tests/ExtendBounds"
	for _, filepath := range GetTestFiles(t, testDir) {
		tc := ReadExtendBoundsTest(filepath)
		t.Run(tc.ID, func(t *testing.T) {
			// Ensure scores slice matches N just in case test file is malformed
			if len(tc.Scores) != tc.N {
				t.Logf("Warning: Test %s Score length (%d) != N (%d)", tc.ID, len(tc.Scores), tc.N)
			}

			l, r := ExtendBounds(tc.Seed, tc.N, tc.Scores, tc.Threshold)

			if l != tc.ExpectedL || r != tc.ExpectedR {
				t.Errorf("Test %s failed.\nSeed: %d, Threshold: %.2f\nGot: (%d, %d)\nExpected: (%d, %d)",
					tc.ID, tc.Seed, tc.Threshold, l, r, tc.ExpectedL, tc.ExpectedR)
			}
		})
	}
}

// --- 4. MarkRegion ---

type MarkRegionTest struct {
	ID        string
	ArrSize   int
	L         int
	R         int
	ExpectedT []int // Indices expected to be true
}

func ReadMarkRegionTest(filepath string) MarkRegionTest {
	data := ReadTestFile(filepath)
	tc := MarkRegionTest{
		ID: data["ID"],
	}
	tc.ArrSize, _ = strconv.Atoi(data["ArrSize"])
	tc.L, _ = strconv.Atoi(data["L"])
	tc.R, _ = strconv.Atoi(data["R"])

	// Parse comma separated indices that should be true
	indicesStr := strings.Split(data["ExpectedTrueIndices"], ",")
	for _, s := range indicesStr {
		val, err := strconv.Atoi(strings.TrimSpace(s))
		if err == nil {
			tc.ExpectedT = append(tc.ExpectedT, val)
		}
	}
	return tc
}

func TestMarkRegion(t *testing.T) {
	testDir := "Tests/MarkRegion"
	for _, filepath := range GetTestFiles(t, testDir) {
		tc := ReadMarkRegionTest(filepath)
		t.Run(tc.ID, func(t *testing.T) {
			arr := make([]bool, tc.ArrSize)
			MarkRegion(arr, tc.L, tc.R)

			// Check consistency
			for i, val := range arr {
				shouldBeTrue := false
				for _, expectedIdx := range tc.ExpectedT {
					if i == expectedIdx {
						shouldBeTrue = true
						break
					}
				}

				if val != shouldBeTrue {
					t.Errorf("Test %s failed at index %d. Got %v, expected %v", tc.ID, i, val, shouldBeTrue)
				}
			}
		})
	}
}

// --- 5. CleanUpStructure ---

type CleanUpTest struct {
	ID       string
	Input    string
	Expected string
}

func ReadCleanUpTest(filepath string) CleanUpTest {
	data := ReadTestFile(filepath)
	tc := CleanUpTest{
		ID:       data["ID"],
		Input:    data["Input"],
		Expected: data["Expected"],
	}
	return tc
}

func TestCleanUpStructure(t *testing.T) {
	testDir := "Tests/CleanUpStructure"
	for _, filepath := range GetTestFiles(t, testDir) {
		tc := ReadCleanUpTest(filepath)
		t.Run(tc.ID, func(t *testing.T) {
			result := CleanUpStructure(tc.Input)
			if result != tc.Expected {
				t.Errorf("Test %s failed.\nInput:    %s\nGot:      %s\nExpected: %s",
					tc.ID, tc.Input, result, tc.Expected)
			}
		})
	}
}

// --- 6. PredictStructure ---
// NOTE: This test relies on Global Maps (hydrophobicity, propensities) being properly defined and exported in your main package.
// If your maps are lowercased (private), this test might fail to compile or run logic correctly unless they are accessible.

type PredictStructureTest struct {
	ID       string
	Sequence string
	Expected string
}

func ReadPredictStructureTest(filepath string) PredictStructureTest {
	data := ReadTestFile(filepath)
	tc := PredictStructureTest{
		ID:       data["ID"],
		Sequence: data["Sequence"],
		Expected: data["Expected"],
	}
	return tc
}

func TestPredictStructure(t *testing.T) {
	testDir := "Tests/PredictStructure"
	for _, filepath := range GetTestFiles(t, testDir) {
		tc := ReadPredictStructureTest(filepath)
		t.Run(tc.ID, func(t *testing.T) {
			// Note: We cannot easily mock the global maps here without changing the source code architecture.
			// Assumes 'hydrophobicity' and 'propensities' are populated in the main package.
			
			result := PredictStructure(tc.Sequence)
			
			// Simple length check first
			if len(result) != len(tc.Expected) {
				t.Fatalf("Test %s length mismatch. Got %d, Expected %d", tc.ID, len(result), len(tc.Expected))
			}

			if result != tc.Expected {
				t.Errorf("Test %s failed.\nSeq: %s\nGot: %s\nExp: %s",
					tc.ID, tc.Sequence, result, tc.Expected)
			}
		})
	}
}

