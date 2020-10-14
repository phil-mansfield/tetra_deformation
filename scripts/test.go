package main

import (
	"os"
	"fmt"
	"math"
	"math/rand"
	"runtime"
	
	"gonum.org/v1/gonum/mat"
	"github.com/phil-mansfield/nbody-utils/box"
	"github.com/phil-mansfield/tetra_deformation/lib"
	"github.com/phil-mansfield/gotetra/render/io"
	"github.com/phil-mansfield/gotetra/render/geom"
	"github.com/phil-mansfield/minnow/go/minh"
)

var (
	//NumCPU = runtime.NumCPU()
	NumCPU = 1

	HlistNames = []string{
		//"../data/hlist_1.00000.minh",
		"",
	}
	SheetNames = []string{
		"../data/sheet000_snap100.dat",
	}
	GridNames = []string{
		"../data/L500_pix250.gtet",
	}
)

const (
	Nbuf = 100
	Npart = 100
)

func initialVolume(L float64, N int) float64 {
	return math.Pow(L/float64(N), 3) / 6
}

func upscaleX(x32 []geom.Vec, Ns, Ng, skip int, x64 [][3]float64) {
	j := 0
	for iz := 0; iz < Ns; iz += skip {
		for iy := 0; iy < Ns; iy += skip {
			for ix := 0; ix < Ns; ix += skip {
				i := ix + iy*Ng + iz*Ng*Ng
				x64[j] = [3]float64{
					float64(x32[i][0]), float64(x32[i][1]), float64(x32[i][2]),
				}
				j++
			}
		}
	}
}

func flagHaloParticles(
	L float64, x [][3]float64, hx [3][]float32, hr []float32, isHalo []bool,
) {
	finder := box.NewFinder(L, x)

	for i := range hr {
		idx := finder.Find(
			[3]float64{float64(hx[0][i]), float64(hx[1][i]), float64(hx[2][i])},
			float64(hr[i]) / 1e3, // kpc -> Mpc
		)
		for _, j := range idx {
			isHalo[j] = true
		}
	}
}

func main() {
	rand.Seed(0)

	runtime.GOMAXPROCS(NumCPU)
	
	order := &lib.Order{1, lib.AllTetra}

	hd := &io.SheetHeader{ }
	io.ReadSheetHeaderAt(SheetNames[0], hd)
	x, tet, def := lib.Buffers(order, SheetNames[0])
	x64 := make([][3]float64, len(def))
	isHalo := make([]bool, len(x64))
	density := make([]float64, len(x64))

	rngBuf := make([][]geom.Vec, Nbuf)
	for i := range rngBuf {
		rngBuf[i] = lib.TetraRngBuf(Npart)
	}
	
	vInit := initialVolume(hd.TotalWidth, 1024 / order.Skip)
	
	for snapIdx := range SheetNames {
		sheet, hlist := SheetNames[snapIdx], HlistNames[snapIdx]
		grid := GridNames[snapIdx]
		
		io.ReadSheetHeaderAt(sheet, hd)
					
		lib.Read(order, sheet, x, tet, def)
		for i := range isHalo {
			isHalo[i] = false
			density[i] = 0
		}

		runtime.GC()
		if hlist != "" {
			fHalo := minh.Open(hlist)
			fCols := fHalo.Floats([]string{"x", "y", "z", "rvir"})
			hx := [3][]float32{fCols["x"], fCols["y"],fCols["z"]}
			hr := fCols["rvir"]

			upscaleX(x, int(hd.SegmentWidth),
				int(hd.GridWidth), order.Skip, x64)
			flagHaloParticles(hd.TotalWidth, x64, hx, hr, isHalo)
			runtime.GC()
		}

		if grid != "" {
			ghd, err := io.ReadGridHeader(grid)
			if err != nil { panic(err.Error()) }
			grid, err := io.ReadGrid(grid)
			if err != nil { panic(err.Error()) }

			lib.GetDensities(ghd, grid, tet, rngBuf, density)
			runtime.GC()
		}
		
		eigBuf, buf := &mat.Eigen{ }, make([]complex128, 3)
		
		w := int(hd.SegmentWidth) / order.Skip
		j := 0
		
		for wz := 0; wz < w; wz++ {
			fname := fmt.Sprintf("../data/grids/grid_L500_%.5f_%d.dat",
				1/(1 + hd.Cosmo.Z), wz)
			f, err := os.Create(fname)
			if err != nil { panic(err.Error()) }
		
			fmt.Fprintln(f, `# 0 - halo flag: 1 -> inside Rvir, 0 -> outside
# 1 - log10(rho / rho_avg)
# 2 - log10(volume / initial volume)
# 3 to 5 - lambda_1 to lambda_3.`)
			for wx := 0; wx < w; wx++ {
				for wy := 0; wy < w; wy++ {
					l1, l2, l3 := lib.DeformationEig(def[j], eigBuf, buf)
					vol := 0.0
					for i := 0; i < 6; i++ {
						vol += tet[i + j*6].Volume()/vInit
					}
				
					flag := 0
					if isHalo[j] { flag = 1 }
					
					fmt.Fprintf(f, "%d %.4f %.4f %.3f %.3f %.3f\n",
						flag, math.Log10(density[j]), math.Log10(vol),
						l1, l2, l3)
					j++
				}
			}
			f.Close()
		}
	}
}
