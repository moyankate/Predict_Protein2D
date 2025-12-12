// Author: Yu-Lun Chen
// Date: 2025-12-11
// Description: Drawing functions for bar charts.
// The drawing functions are mainly derived from AI tools.
// Comments are mainly created by AI tools

package main

import (
    "canvas"
    "image/color"
    "math"
)


// ----- Draw Barchart ----- //

// getColorRGBA Convert structure character to an RGBA value acceptable by canvas.MakeColor
func getColorRGBA(structureChar rune) color.Color {
    switch structureChar {
    case 'H':
        // Red (Alpha Helix)
        return canvas.MakeColor(255, 0, 0)
    case 'E':
        // Yellow (Beta Sheet)
        return canvas.MakeColor(255, 215, 0)
    case 'C':
        // Blue (Coil / Unpredicted)
        return canvas.MakeColor(25, 25, 112)
    default:
        // Grey (Fallback for unexpected chars)
        return canvas.MakeColor(128, 128, 128)
    }
}

// drawPredictionBarChartCanvas draws the Chou-Fasman prediction result as a PNG image
// The input is the final predicted structure sequence (H/E/C).
func drawPredictionBarChartCanvas(predictedSequence string, filename string) error {
    const (
        barWidth    = 3        // Width of each residue bar
        barHeight   = 20       // Height of each residue bar
        padding     = 20       // Padding around the edges
        lineHeight  = 25       // Height between lines
        charsPerLine = 100     // Number of residues displayed per line
        tickSize    = 1        // Width of the simulated tick
        tickLength  = 5        // Length of the simulated tick
    )

    seqLen := len(predictedSequence) // Use the length of the predicted sequence
    numLines := int(math.Ceil(float64(seqLen) / float64(charsPerLine)))
    
    imgWidth := padding*2 + charsPerLine*barWidth
    imgHeight := padding*2 + numLines*(lineHeight + barHeight) + 100 

    // 1. Create canvas
    c := canvas.CreateNewCanvas(imgWidth, imgHeight)

    // Set background to white
    c.SetFillColor(canvas.MakeColor(255, 255, 255))
    c.Clear()
    
    // Set default drawing color to black
    c.SetStrokeColor(canvas.MakeColor(0, 0, 0)) 
    c.SetFillColor(canvas.MakeColor(0, 0, 0)) 
    
    // 2. Iterate and draw bars for each residue
    for i, structureChar := range predictedSequence {
        lineNum := i / charsPerLine
        charInLine := i % charsPerLine

        x := float64(padding + charInLine*barWidth)
        y := float64(padding + 50 + lineNum*(lineHeight + barHeight)) 

        barColor := getColorRGBA(structureChar) // Get color based on the structure character
        
        // Draw the main structure bar (rectangle)
        c.SetFillColor(barColor)
        
        // Define rectangle path
        c.MoveTo(x, y)
        c.LineTo(x + barWidth, y)
        c.LineTo(x + barWidth, y + barHeight)
        c.LineTo(x, y + barHeight)
        
        // Fill the rectangle
        c.Fill() 
        
        // Draw simulated tick marks (every 10 residues)
        if i % 10 == 0 {
            // Define a small vertical rectangle as the tick mark
            tickX := x
            tickY := y + barHeight
            
            c.SetFillColor(canvas.MakeColor(0, 0, 0)) // Black for tick
            
            // Define tick mark path
            c.MoveTo(tickX, tickY)
            c.LineTo(tickX + tickSize, tickY)
            c.LineTo(tickX + tickSize, tickY + tickLength)
            c.LineTo(tickX, tickY + tickLength)
            
            c.Fill() 
        }
    }

    // 3. Save as PNG file
    c.SaveToPNG(filename) 
    
    return nil
}