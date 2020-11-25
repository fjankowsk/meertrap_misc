package main

import (
	"bufio"
	"encoding/csv"
	"fmt"
	"log"
	"os"
	//"sort"
	"strconv"
)

type Data struct {
	x float64
	y float64
	dist float64
}

func load_data(filename string) ([][]float64, error) {
	f, err := os.Open(filename)

	if err != nil {
		error := fmt.Errorf("Could not open file: %s, %s", filename, err)
		return nil, error
	}

	defer f.Close()

	reader := csv.NewReader(bufio.NewReader(f))
	reader.Comma = '\t'

	lines, err := reader.ReadAll()

	if err != nil {
		error := fmt.Errorf("Could not parse csv data: %s", err)
		return nil, error
	}

	var data [][]float64

	for _, line := range lines {
		x, _ := strconv.ParseFloat(line[0], 64)
		y, _ := strconv.ParseFloat(line[1], 64)

		item := []float64{x, y}

		data = append(data, item)
	}

	return data, nil
}


//func get_beam_packing(data [][]float64){
//	const nbeams = 396
//	const bunch = 6
//
//	sort.Float64s(data)
//}


func main() {
	infile := "input/134.0696_0.0_beam_pos.dat"

	data, err := load_data(infile)
	if err != nil {
		log.Fatalf("Could not load data from file: %s, %s", infile, err)
	}

	for _, blo := range data {
		fmt.Println(blo[0], blo[1])
	}
}