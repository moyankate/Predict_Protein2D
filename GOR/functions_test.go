// Author: Sorro Sun
// Date: 2025-12-11
// Description: Unit tests for GOR secondary structure prediction functions in Go
// Note: Test cases are stored in separate .txt files under test/ directory.

package gor

import (
	"bufio"
	"encoding/json"
	"fmt"
	"math"
	"os"
	"sort"
	"strconv"
	"strings"
	"testing"
)

// --- General Helper Functions ---

// ReadTestFile loads the file content into a map for structured parsing.
func ReadTestFile(filepath string) map[string]string {
	f, err := os.Open(filepath)
	if err != nil {
		panic(fmt.Sprintf("Error reading file %s: %v", filepath, err))
	}
	defer f.Close()

	data := make(map[string]string)
	scanner := bufio.NewScanner(f)

	// Custom parsing to handle multi-line values (like matrices/models in JSON)
	var currentKey string
	var currentValBuilder strings.Builder

	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())

		// If line starts with "#", it's a new key
		if strings.HasPrefix(line, "#") {
			// Save previous key if exists
			if currentKey != "" {
				data[currentKey] = strings.TrimSpace(currentValBuilder.String())
				currentValBuilder.Reset()
			}

			// Parse new key
			parts := strings.SplitN(strings.TrimPrefix(line, "# "), ":", 2)
			if len(parts) >= 1 {
				currentKey = strings.TrimSpace(parts[0])
				if len(parts) == 2 && strings.TrimSpace(parts[1]) != "" {
					// Inline value
					currentValBuilder.WriteString(strings.TrimSpace(parts[1]))
				}
			}
		} else {
			// Continuation of the current key's value
			if currentKey != "" {
				currentValBuilder.WriteString(line)
				currentValBuilder.WriteString("\n")
			}
		}
	}
	// Save last key
	if currentKey != "" {
		data[currentKey] = strings.TrimSpace(currentValBuilder.String())
	}

	return data
}

// GetTestFiles reads all .txt files in the given directory.
func GetTestFiles(t *testing.T, dir string) []string {
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

// createTempFile creates a file with content for testing IO functions
func createTempFile(t *testing.T, content string) string {
	tmpFile, err := os.CreateTemp("", "gortest_*.txt")
	if err != nil {
		t.Fatalf("Failed to create temp file: %v", err)
	}
	if _, err := tmpFile.WriteString(content); err != nil {
		t.Fatalf("Failed to write to temp file: %v", err)
	}
	tmpFile.Close()
	return tmpFile.Name()
}

// --- 1. ParsePSSM ---

type ParsePSSMTest struct {
	ID          string
	FileContent string
	ExpectedLen int
	CheckValues string // e.g., "0:0.01, 1:0.45" to check specific positions
	ShouldErr   bool
}

func ReadParsePSSMTest(filepath string) ParsePSSMTest {
	data := ReadTestFile(filepath)
	tc := ParsePSSMTest{
		ID:          data["ID"],
		FileContent: data["FileContent"],
		CheckValues: data["CheckValues"],
		ShouldErr:   data["ShouldErr"] == "true",
	}
	tc.ExpectedLen, _ = strconv.Atoi(data["ExpectedLen"])
	return tc
}

func TestParsePSSM(t *testing.T) {
	testDir := "test/ParsePSSM"
	for _, filepath := range GetTestFiles(t, testDir) {
		tc := ReadParsePSSMTest(filepath)
		t.Run(tc.ID, func(t *testing.T) {
			tmpPath := createTempFile(t, tc.FileContent)
			defer os.Remove(tmpPath)

			result, err := ParsePSSM(tmpPath)

			if tc.ShouldErr {
				if err == nil {
					t.Errorf("Expected error but got none")
				}
			} else {
				if err != nil {
					t.Errorf("Unexpected error: %v", err)
				} else {
					if len(result) != tc.ExpectedLen {
						t.Errorf("Length mismatch. Got %d, Expected %d", len(result), tc.ExpectedLen)
					}
					// Check specific values if provided
					if tc.CheckValues != "" {
						pairs := strings.Split(tc.CheckValues, ",")
						for _, pair := range pairs {
							parts := strings.Split(strings.TrimSpace(pair), ":")
							if len(parts) == 2 {
								rowIdx, _ := strconv.Atoi(parts[0])
								expectedVal, _ := strconv.ParseFloat(parts[1], 64)
								if rowIdx < len(result) && len(result[rowIdx]) > 0 {
									// Check first element for simplicity
									if math.Abs(result[rowIdx][0]-expectedVal) > 1e-6 {
										t.Errorf("Row %d: Got %f, Expected %f", rowIdx, result[rowIdx][0], expectedVal)
									}
								}
							}
						}
					}
				}
			}
		})
	}
}

// --- 2. ParseDSSP ---

type ParseDSSPTest struct {
	ID          string
	FileContent string
	Expected    string
	ShouldErr   bool
}

func ReadParseDSSPTest(filepath string) ParseDSSPTest {
	data := ReadTestFile(filepath)
	return ParseDSSPTest{
		ID:          data["ID"],
		FileContent: data["FileContent"],
		Expected:    data["Expected"],
		ShouldErr:   data["ShouldErr"] == "true",
	}
}

func TestParseDSSP(t *testing.T) {
	testDir := "test/ParseDSSP"
	for _, filepath := range GetTestFiles(t, testDir) {
		tc := ReadParseDSSPTest(filepath)
		t.Run(tc.ID, func(t *testing.T) {
			tmpPath := createTempFile(t, tc.FileContent)
			defer os.Remove(tmpPath)

			result, err := ParseDSSP(tmpPath)

			if tc.ShouldErr {
				if err == nil {
					t.Errorf("Expected error but got none")
				}
			} else {
				if err != nil {
					t.Errorf("Unexpected error: %v", err)
				}
				if result != tc.Expected {
					t.Errorf("Got: %q, Expected: %q", result, tc.Expected)
				}
			}
		})
	}
}

// --- 3. ParseFASTA ---

type ParseFASTATest struct {
	ID          string
	FileContent string
	Expected    string
	ShouldErr   bool
}

func ReadParseFASTATest(filepath string) ParseFASTATest {
	data := ReadTestFile(filepath)
	return ParseFASTATest{
		ID:          data["ID"],
		FileContent: data["FileContent"],
		Expected:    data["Expected"],
		ShouldErr:   data["ShouldErr"] == "true",
	}
}

func TestParseFASTA(t *testing.T) {
	testDir := "test/ParseFASTA"
	for _, filepath := range GetTestFiles(t, testDir) {
		tc := ReadParseFASTATest(filepath)
		t.Run(tc.ID, func(t *testing.T) {
			tmpPath := createTempFile(t, tc.FileContent)
			defer os.Remove(tmpPath)

			result, err := ParseFASTA(tmpPath)

			if tc.ShouldErr {
				if err == nil {
					t.Errorf("Expected error but got none")
				}
			} else {
				if err != nil {
					t.Errorf("Unexpected error: %v", err)
				}
				if result != tc.Expected {
					t.Errorf("Got: %s, Expected: %s", result, tc.Expected)
				}
			}
		})
	}
}

// --- 4. SeqToProfile ---

type SeqToProfileTest struct {
	ID             string
	Sequence       string
	ExpectedLen    int
	ExpectedValues map[int]string // Map index -> expected AA
}

func ReadSeqToProfileTest(filepath string) SeqToProfileTest {
	data := ReadTestFile(filepath)

	tc := SeqToProfileTest{
		ID:          data["ID"],
		Sequence:    data["Sequence"],
		ExpectedLen: 0,
	}
	tc.ExpectedLen, _ = strconv.Atoi(data["ExpectedLen"])

	tc.ExpectedValues = make(map[int]string)
	if vals, ok := data["CheckValues"]; ok && vals != "" {
		parts := strings.Split(vals, ",")
		for _, p := range parts {
			kv := strings.Split(strings.TrimSpace(p), ":")
			if len(kv) == 2 {
				idx, _ := strconv.Atoi(kv[0])
				tc.ExpectedValues[idx] = kv[1]
			}
		}
	}
	return tc
}

func TestSeqToProfile(t *testing.T) {
	aaList := "ARNDCQEGHILKMFPSTWYV"

	testDir := "test/SeqToProfile"
	for _, filepath := range GetTestFiles(t, testDir) {
		tc := ReadSeqToProfileTest(filepath)
		t.Run(tc.ID, func(t *testing.T) {
			profile := SeqToProfile(tc.Sequence)

			if len(profile) != tc.ExpectedLen {
				t.Errorf("Length mismatch. Got %d, Expected %d", len(profile), tc.ExpectedLen)
			}

			for rowIdx, expectedAA := range tc.ExpectedValues {
				if rowIdx >= len(profile) {
					continue
				}
				found := -1
				for colIdx, val := range profile[rowIdx] {
					if math.Abs(val-1.0) < 1e-9 {
						found = colIdx
						break
					}
				}

				expectedIndex := strings.IndexByte(aaList, expectedAA[0])
				if found != expectedIndex {
					t.Errorf("Row %d: Expected active bit for %s (index %d), but found active at index %d",
						rowIdx, expectedAA, expectedIndex, found)
				}
			}
		})
	}
}

// --- 5. sumProfile ---

type SumProfileTest struct {
	ID          string
	ProfileJSON string
	Expected    float64
}

func ReadSumProfileTest(filepath string) SumProfileTest {
	data := ReadTestFile(filepath)
	tc := SumProfileTest{
		ID:       data["ID"],
		Expected: 0.0,
	}
	tc.ProfileJSON = data["Profile"]
	tc.Expected, _ = strconv.ParseFloat(data["Expected"], 64)
	return tc
}

func TestSumProfile(t *testing.T) {
	testDir := "test/sumProfile"
	for _, filepath := range GetTestFiles(t, testDir) {
		tc := ReadSumProfileTest(filepath)
		t.Run(tc.ID, func(t *testing.T) {
			var profile [][]float64
			if err := json.Unmarshal([]byte(tc.ProfileJSON), &profile); err != nil {
				t.Fatalf("Invalid Profile JSON in test file: %v", err)
			}

			result := sumProfile(profile)

			if math.Abs(result-tc.Expected) > 1e-6 {
				t.Errorf("Got: %f, Expected: %f", result, tc.Expected)
			}
		})
	}
}

// --- 6. TrainGOR ---

type TrainGORTest struct {
	ID         string
	IDList     string
	PSSMFile   string // PSSM file content for each ID
	DSSPFile   string // DSSP file content for each ID
	WindowSize int
	ShouldErr  bool
}

func ReadTrainGORTest(filepath string) TrainGORTest {
	data := ReadTestFile(filepath)
	tc := TrainGORTest{
		ID:         data["ID"],
		IDList:     data["IDList"],
		PSSMFile:   data["PSSMFile"],
		DSSPFile:   data["DSSPFile"],
		WindowSize: 3,
		ShouldErr:  data["ShouldErr"] == "true",
	}
	if ws, err := strconv.Atoi(data["WindowSize"]); err == nil {
		tc.WindowSize = ws
	}
	return tc
}

func TestTrainGOR(t *testing.T) {
	testDir := "test/TrainGOR"
	for _, filepath := range GetTestFiles(t, testDir) {
		tc := ReadTrainGORTest(filepath)
		t.Run(tc.ID, func(t *testing.T) {
			// Create temp directories
			pssmDir, err := os.MkdirTemp("", "pssm_*")
			if err != nil {
				t.Fatalf("Failed to create temp PSSM dir: %v", err)
			}
			defer os.RemoveAll(pssmDir)

			dsspDir, err := os.MkdirTemp("", "dssp_*")
			if err != nil {
				t.Fatalf("Failed to create temp DSSP dir: %v", err)
			}
			defer os.RemoveAll(dsspDir)

			// Parse ID list and create files
			ids := strings.Split(strings.TrimSpace(tc.IDList), "\n")
			var validIDs []string
			for _, id := range ids {
				id = strings.TrimSpace(id)
				if id == "" {
					continue
				}
				validIDs = append(validIDs, id)

				// Create PSSM file
				if tc.PSSMFile != "" {
					pssmPath := fmt.Sprintf("%s/%s.pssm", pssmDir, id)
					pssmF, err := os.Create(pssmPath)
					if err == nil {
						pssmF.WriteString(tc.PSSMFile)
						pssmF.Close()
					}
				}

				// Create DSSP file
				if tc.DSSPFile != "" {
					dsspPath := fmt.Sprintf("%s/%s.dssp", dsspDir, id)
					dsspF, err := os.Create(dsspPath)
					if err == nil {
						dsspF.WriteString(tc.DSSPFile)
						dsspF.Close()
					}
				}
			}

			// Create ID list file
			idListPath := createTempFile(t, tc.IDList)
			defer os.Remove(idListPath)

			model, err := TrainGOR(idListPath, pssmDir, dsspDir, tc.WindowSize)

			if tc.ShouldErr {
				if err == nil {
					t.Errorf("Expected error but got none")
				}
			} else {
				if err != nil {
					t.Errorf("Unexpected error: %v", err)
				} else if model == nil {
					t.Errorf("Got nil model")
				} else {
					if model.WindowSize != tc.WindowSize {
						t.Errorf("WindowSize mismatch. Got %d, Expected %d", model.WindowSize, tc.WindowSize)
					}
				}
			}
		})
	}
}

// --- 7. SaveGORModel ---

type SaveGORModelTest struct {
	ID        string
	ModelJSON string
	ShouldErr bool
}

func ReadSaveGORModelTest(filepath string) SaveGORModelTest {
	data := ReadTestFile(filepath)
	return SaveGORModelTest{
		ID:        data["ID"],
		ModelJSON: data["Model"],
		ShouldErr: data["ShouldErr"] == "true",
	}
}

func TestSaveGORModel(t *testing.T) {
	testDir := "test/SaveGORModel"
	for _, filepath := range GetTestFiles(t, testDir) {
		tc := ReadSaveGORModelTest(filepath)
		t.Run(tc.ID, func(t *testing.T) {
			var model GORModel
			if err := json.Unmarshal([]byte(tc.ModelJSON), &model); err != nil {
				t.Fatalf("Invalid Model JSON in test file: %v", err)
			}

			tmpPath := createTempFile(t, "")
			os.Remove(tmpPath) // Remove empty file
			tmpPath = tmpPath + ".json"

			err := SaveGORModel(&model, tmpPath)
			defer os.Remove(tmpPath)

			if tc.ShouldErr {
				if err == nil {
					t.Errorf("Expected error but got none")
				}
			} else {
				if err != nil {
					t.Errorf("Unexpected error: %v", err)
				} else {
					// Verify file exists and can be loaded
					loaded, err := LoadGORModel(tmpPath)
					if err != nil {
						t.Errorf("Failed to load saved model: %v", err)
					} else if loaded.WindowSize != model.WindowSize {
						t.Errorf("WindowSize mismatch after save/load")
					}
				}
			}
		})
	}
}

// --- 8. LoadGORModel ---

type LoadGORModelTest struct {
	ID         string
	ModelFile  string
	ExpectedWS int
	ShouldErr  bool
}

func ReadLoadGORModelTest(filepath string) LoadGORModelTest {
	data := ReadTestFile(filepath)
	tc := LoadGORModelTest{
		ID:         data["ID"],
		ModelFile:  data["ModelFile"],
		ExpectedWS: 3,
		ShouldErr:  data["ShouldErr"] == "true",
	}
	if ws, err := strconv.Atoi(data["ExpectedWS"]); err == nil {
		tc.ExpectedWS = ws
	}
	return tc
}

func TestLoadGORModel(t *testing.T) {
	testDir := "test/LoadGORModel"
	for _, filepath := range GetTestFiles(t, testDir) {
		tc := ReadLoadGORModelTest(filepath)
		t.Run(tc.ID, func(t *testing.T) {
			// Create temp file from ModelFile content
			tmpPath := createTempFile(t, tc.ModelFile)
			defer os.Remove(tmpPath)

			result, err := LoadGORModel(tmpPath)

			if tc.ShouldErr {
				if err == nil {
					t.Errorf("Expected error but got none")
				}
			} else {
				if err != nil {
					t.Errorf("Unexpected error: %v", err)
				} else if result == nil {
					t.Errorf("Got nil model")
				} else {
					if result.WindowSize != tc.ExpectedWS {
						t.Errorf("WindowSize mismatch. Got %d, Expected %d", result.WindowSize, tc.ExpectedWS)
					}
				}
			}
		})
	}
}

// --- 9. PredictGOR ---

type PredictGORTest struct {
	ID          string
	ModelJSON   string
	ProfileJSON string
	Expected    string
}

func ReadPredictGORTest(filepath string) PredictGORTest {
	data := ReadTestFile(filepath)
	return PredictGORTest{
		ID:          data["ID"],
		ModelJSON:   data["Model"],
		ProfileJSON: data["Profile"],
		Expected:    data["Expected"],
	}
}

func TestPredictGOR(t *testing.T) {
	testDir := "test/PredictGOR"
	for _, filepath := range GetTestFiles(t, testDir) {
		tc := ReadPredictGORTest(filepath)
		t.Run(tc.ID, func(t *testing.T) {
			var model GORModel
			if err := json.Unmarshal([]byte(tc.ModelJSON), &model); err != nil {
				t.Fatalf("Invalid Model JSON in test file: %v", err)
			}

			var profile [][]float64
			if err := json.Unmarshal([]byte(tc.ProfileJSON), &profile); err != nil {
				t.Fatalf("Invalid Profile JSON in test file: %v", err)
			}

			result, err := PredictGOR(&model, profile)
			if err != nil {
				t.Errorf("Unexpected error: %v", err)
			}

			if result != tc.Expected {
				t.Errorf("Got: %s, Expected: %s", result, tc.Expected)
			}
		})
	}
}
