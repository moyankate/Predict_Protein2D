// Author: Yu-Lun Chen
// Date: 2025-12-11
// Description: Definition of datatypes using in the protein secondary prediction project.

package main

// Global variable
var SEQUENCE string

// Region: to store the start and end position for resulting regions
type Region struct {
	Start int
	End int
}
