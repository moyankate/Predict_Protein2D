package gor

import (
	"bufio"
	"encoding/json"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"strconv"
	"strings"
)

// GORModel holds the trained probability tables.
type GORModel struct {
	WindowSize int         `json:"window_size"`
	AAList     []string    `json:"aa_list"`
	H          [][]float64 `json:"H"` // helix
	E          [][]float64 `json:"E"` // strand
	C          [][]float64 `json:"C"` // coil ("-")
}

//
// ===== util-like helpers (Python utils.py) =====
//

// ParsePSSM reads a PSSM file and returns a profile [L][20].
func ParsePSSM(filename string) ([][]float64, error) {
	f, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	defer f.Close()

	var lines []string
	scanner := bufio.NewScanner(f)
	for scanner.Scan() {
		lines = append(lines, scanner.Text())
	}
	if err := scanner.Err(); err != nil {
		return nil, err
	}

	if len(lines) <= 9 {
		return [][]float64{}, nil
	}
	// mimic Python lines[3:-6]
	lines = lines[3 : len(lines)-6]

	var profile [][]float64
	for _, line := range lines {
		fields := strings.Fields(line)
		if len(fields) < 24 {
			continue
		}
		cols := fields[22 : len(fields)-2] // [22:-2]
		row := make([]float64, len(cols))
		for i, s := range cols {
			v, err := strconv.ParseFloat(s, 64)
			if err != nil {
				return nil, fmt.Errorf("ParsePSSM: %v", err)
			}
			row[i] = v / 100.0
		}
		profile = append(profile, row)
	}
	return profile, nil
}

// ParseDSSP reads the DSSP SS string (second line).
func ParseDSSP(filename string) (string, error) {
	f, err := os.Open(filename)
	if err != nil {
		return "", err
	}
	defer f.Close()

	r := bufio.NewReader(f)
	// first line
	if _, err := r.ReadString('\n'); err != nil {
		return "", err
	}
	// second line
	line, err := r.ReadString('\n')
	if err != nil && err != io.EOF {
		return "", err
	}
	return strings.TrimSpace(line), nil
}

// ParseFASTA reads a simple FASTA (header + one sequence line).
func ParseFASTA(filename string) (string, error) {
	f, err := os.Open(filename)
	if err != nil {
		return "", err
	}
	defer f.Close()

	scanner := bufio.NewScanner(f)
	// header
	if !scanner.Scan() {
		return "", fmt.Errorf("empty FASTA: %s", filename)
	}
	// sequence
	if !scanner.Scan() {
		return "", fmt.Errorf("no sequence line in FASTA: %s", filename)
	}
	return strings.TrimSpace(scanner.Text()), nil
}

// SeqToProfile converts sequence -> one-hot [L][20].
func SeqToProfile(seq string) [][]float64 {
	aaList := []rune{'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I',
		'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'}

	index := make(map[rune]int)
	for i, aa := range aaList {
		index[aa] = i
	}

	profile := make([][]float64, len(seq))
	for i, r := range seq {
		row := make([]float64, len(aaList))
		if j, ok := index[r]; ok {
			row[j] = 1.0
		}
		profile[i] = row
	}
	return profile
}

func sumProfile(profile [][]float64) float64 {
	var s float64
	for _, row := range profile {
		for _, v := range row {
			s += v
		}
	}
	return s
}

//
// ===== training =====
//

// TrainGOR trains a GOR model from:
//  - idListFile: file with protein IDs
//  - pssmDir: directory with <id>.pssm
//  - dsspDir: directory with <id>.dssp
func TrainGOR(idListFile, pssmDir, dsspDir string, windowSize int) (*GORModel, error) {
	aaList := []string{"A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
		"L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"}

	// allocate matrices [windowSize][20]
	newMat := func() [][]float64 {
		m := make([][]float64, windowSize)
		for i := range m {
			m[i] = make([]float64, len(aaList))
		}
		return m
	}

	model := &GORModel{
		WindowSize: windowSize,
		AAList:     aaList,
		H:          newMat(),
		E:          newMat(),
		C:          newMat(),
	}

	f, err := os.Open(idListFile)
	if err != nil {
		return nil, err
	}
	defer f.Close()
	scanner := bufio.NewScanner(f)

	halfWindow := (windowSize - 1) / 2

	for scanner.Scan() {
		id := strings.TrimSpace(scanner.Text())
		if id == "" {
			continue
		}

		pssmPath := filepath.Join(pssmDir, id+".pssm")
		dsspPath := filepath.Join(dsspDir, id+".dssp")

		profile, err := ParsePSSM(pssmPath)
		if err != nil {
			return nil, fmt.Errorf("ParsePSSM %s: %w", pssmPath, err)
		}
		if len(profile) == 0 || sumProfile(profile) == 0 {
			continue
		}

		ss, err := ParseDSSP(dsspPath)
		if err != nil {
			return nil, fmt.Errorf("ParseDSSP %s: %w", dsspPath, err)
		}
		// replace C -> -
		ss = strings.ReplaceAll(ss, "C", "-")
		if len(ss) != len(profile) {
			return nil, fmt.Errorf("length mismatch %s: profile=%d ss=%d",
				id, len(profile), len(ss))
		}

		for i := 0; i < len(profile); i++ {
			var mat [][]float64
			switch ss[i] {
			case 'H':
				mat = model.H
			case 'E':
				mat = model.E
			default:
				mat = model.C
			}

			start := i - halfWindow
			if start < 0 {
				start = 0
			}
			end := i + halfWindow + 1
			if end > len(profile) {
				end = len(profile)
			}

			for j := start; j < end; j++ {
				rel := j - i + halfWindow // 0..windowSize-1
				for k := 0; k < len(profile[i]); k++ {
					mat[rel][k] += profile[j][k]
				}
			}
		}
	}
	if err := scanner.Err(); err != nil {
		return nil, err
	}

	normalize := func(m [][]float64) {
		for r := range m {
			var sum float64
			for _, v := range m[r] {
				sum += v
			}
			if sum > 0 {
				for c := range m[r] {
					m[r][c] /= sum
				}
			}
		}
	}

	normalize(model.H)
	normalize(model.E)
	normalize(model.C)

	return model, nil
}

// SaveGORModel saves model as JSON.
func SaveGORModel(model *GORModel, filename string) error {
	f, err := os.Create(filename)
	if err != nil {
		return err
	}
	defer f.Close()

	enc := json.NewEncoder(f)
	enc.SetIndent("", "  ")
	return enc.Encode(model)
}

// LoadGORModel loads model from JSON.
func LoadGORModel(filename string) (*GORModel, error) {
	f, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	defer f.Close()

	var m GORModel
	if err := json.NewDecoder(f).Decode(&m); err != nil {
		return nil, err
	}
	return &m, nil
}

//
// ===== prediction =====
//

// PredictGOR predicts SS string ("HE-..") from a profile [L][20].
func PredictGOR(model *GORModel, profile [][]float64) (string, error) {
	if len(profile) == 0 {
		return "", nil
	}
	windowSize := model.WindowSize
	halfWindow := (windowSize - 1) / 2

	out := make([]rune, len(profile))

	for i := 0; i < len(profile); i++ {
		var scoreH, scoreE, scoreC float64

		start := i - halfWindow
		if start < 0 {
			start = 0
		}
		end := i + halfWindow + 1
		if end > len(profile) {
			end = len(profile)
		}

		for j := start; j < end; j++ {
			rel := j - i + halfWindow
			for k := 0; k < len(profile[i]); k++ {
				scoreH += profile[j][k] * model.H[rel][k]
				scoreE += profile[j][k] * model.E[rel][k]
				scoreC += profile[j][k] * model.C[rel][k]
			}
		}

		label := 'H'
		maxScore := scoreH
		if scoreE > maxScore {
			maxScore = scoreE
			label = 'E'
		}
		if scoreC > maxScore {
			label = 'C'
		}
		out[i] = label
	}

	return string(out), nil
}
