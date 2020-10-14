package main

import (
	"encoding/binary"
	"fmt"
	"math"
	"os"
	"runtime"

	"gonum.org/v1/gonum/mat"
	
	"github.com/phil-mansfield/gotetra/render/geom"
	"github.com/phil-mansfield/gotetra/render/io"
	
	"github.com/phil-mansfield/minnow/go/minh"

	"github.com/phil-mansfield/tetra_deformation/lib"
)

const (
	NumCPU = 1
	SnapMin, SnapMax = 98, 100
	HaloDir = "/home/phil/code/src/github.com/phil-mansfield/tetra_deformation/data/L63/haloes"
	ParticleDir = "/home/phil/code/src/github.com/phil-mansfield/tetra_deformation/data/L63/particles"
	InfoDir = "/home/phil/code/src/github.com/phil-mansfield/tetra_deformation/data/L63/Lagrangian_info"
)

var (
	CellMin, CellMax = []int{ 0, 0, 0 }, []int{ 7, 7, 7 }
)

type Halo struct {
	Desc *Halo
	Prog []*Halo
	Idx int
}

type IDs struct {
	N int
	ID, DescID []int64
}

func ReadIDs(file string) *IDs {
	runtime.GC()
	f := minh.Open(file)
	defer f.Close()
	
	iCols := f.Ints([]string{ "id", "desc_id" })
	
	return &IDs{
		N: len(iCols["id"]),
		ID: iCols["id"],
		DescID: iCols["desc_id"],
	}
}


func UnlinkedHaloes(h *IDs) []Halo {
	haloes := make([]Halo, h.N)
	for i := range haloes {
		haloes[i].Idx = i
	}

	return haloes
}

func LinkHaloes(pH, dH []Halo, pID, dID *IDs) {
	// Find ID range so we can construct a lookup table.
	min, max := dID.ID[0], dID.ID[0]
	for _, id := range dID.ID {
		if id < min {
			min = id
		} else if id > max {
			max = id
		}
	}
	
	// Make lookup table.
	lookup := make([]int, max - min + 1)
	for d := range dID.ID {
		lookup[dID.ID[d] - min] = d
	}

	n := make([]int, len(lookup))
	
	// Check how many mergers there are for each object so we can avoid
	// heap thrashing with append calls.
	for p := range pID.ID {
		if pID.DescID[p] < 0 { continue }
		d := lookup[pID.DescID[p] - min]
		n[d]++
	}

	for d := range dH {
		if n[d] > 0 {
			dH[d].Prog = make([]*Halo, 0, n[d])
		}
	}


	// Link objects using the lookup table.
	for p := range pID.ID {
		if pID.DescID[p] < 0 { continue }
		
		d := lookup[pID.DescID[p] - min]
		pH[p].Desc = &dH[d]
		dH[d].Prog = append(dH[d].Prog, &pH[p])
	}
}


func ReadMergerTree() [][]Halo {
	pID, dID := &IDs{ }, &IDs{ }

	haloes := make([][]Halo, SnapMax + 1)
	hd := &io.SheetHeader{ }
	
	for snap := SnapMax; snap >= SnapMin; snap-- {
		runtime.GC()
		
		pFile := fmt.Sprintf("%s/snapdir_%03d/sheet000.dat", HaloDir, snap)
		io.ReadSheetHeaderAt(pFile, hd)
		scale := 1/(hd.Cosmo.Z + 1)
		hFile := fmt.Sprintf("%s/hlist_%.5f.minh", HaloDir, scale)
		
		pID = ReadIDs(hFile)
		pH := UnlinkedHaloes(pID)
		
		haloes[snap] = pH
		if snap != SnapMax {
			LinkHaloes(haloes[snap], haloes[snap-1], pID, dID)
		}
		dID = pID
	}

	return haloes
}

func LagrangianDensity(
	cx, cy, cz int,
	order *lib.Order, x []geom.Vec, tet []geom.Tetra, def []*mat.Dense,
	density []float32, lambda [3][]float32,
) {
	hd := &io.SheetHeader{ }
	icFile := fmt.Sprintf("%s/IC/sheet%d%d%d.dat", ParticleDir, cx, cy, cz)
	
	io.ReadSheetHeaderAt(icFile, hd)
	lib.Read(order, icFile, x, tet, def)

	initVol := math.Pow(hd.TotalWidth / float64(hd.CountWidth), 3)
	N := len(def)
	buf, eigen := []complex128{ 0, 0, 0 }, &mat.Eigen{ }
	
	for i := 0; i < N; i++ {
		vol := 0.0
		for j := 0; j < 6; j++ {
			vol += tet[i*6 + j].Volume()
		}

		density[i] = float32(initVol / vol)
		l1, l2, l3 := lib.DeformationEig(def[i], eigen, buf)
		lambda[0][i] = float32(l1)
		lambda[1][i] = float32(l2)
		lambda[2][i] = float32(l3)
	}
	
	outFile := fmt.Sprintf("%s/L_density%d%d%d.dat", InfoDir, cx, cy, cz)
	f, err := os.Create(outFile)
	if err != nil { err.Error() }
	defer f.Close()

	binary.Write(f, binary.LittleEndian, density)
	binary.Write(f, binary.LittleEndian, lambda[0])
	binary.Write(f, binary.LittleEndian, lambda[1])
	binary.Write(f, binary.LittleEndian, lambda[2])
}

func FirstAccretion(cx, cy, cz int, tree [][]Halo, x []geom.Vec) {
}

func FinalRadius(cx, cy, cz int, x []geom.Vec) {
}

func main() {
	runtime.GOMAXPROCS(NumCPU)

	order := &lib.Order{ 1, lib.AllTetra }
	defaultFile := fmt.Sprintf(
		"%s/snapdir_%03d/sheet000.dat", ParticleDir, SnapMax,
	)
	x, tet, def := lib.Buffers(order, defaultFile)
	tree := ReadMergerTree()

	N := len(def)
	density := make([]float32, N)
	lambda := [3][]float32{
		make([]float32, N), make([]float32, N), make([]float32, N),
	}
	
	for cx := CellMin[0]; cx <= CellMax[0]; cx++ {
		for cy := CellMin[1]; cy <= CellMax[1]; cy++ {
			for cz := CellMin[2]; cz <= CellMax[2]; cz++ {
				LagrangianDensity(
					cx, cy, cz, order,
					x, tet, def,
					density, lambda,
				)
				FirstAccretion(cx, cy, cz, tree, x)
				FinalRadius(cx, cy, cz, x)
			}
		}
	}
}
