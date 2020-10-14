package lib

import (
	"fmt"
	"math"
	"math/rand"
	"runtime"

	"gonum.org/v1/gonum/mat"
	
	gtet_rand "github.com/phil-mansfield/gotetra/math/rand"
	"github.com/phil-mansfield/gotetra/render/geom"
	"github.com/phil-mansfield/gotetra/render/io"
	
	"github.com/phil-mansfield/nbody-utils/thread"
)

// GetDensities writes densities associated with an array of tetrahedra into the
// array out. These are interpolated off the grid using CIC. Requires buffers of
// randompoint sampled from within a tetrahedra. (i.e. the output of
// TetraRngBuf()) This function is multithreaded and will automatically run
// the hot loop on GOMAXPROCS cores.
func GetDensities(
	ghd *io.GridHeader, grid[]float64, tet []geom.Tetra,
	rngBuf [][]geom.Vec, out []float64,
) {
	numCPU := runtime.GOMAXPROCS(-1)
	thread.SplitArray(len(out), numCPU, func(worker, start, end, step int)  {
		vecBuf := make([]geom.Vec, len(rngBuf[0]))
		for i := start; i < end; i += step {
			for j := 0; j < 6; j++ {
				bufIdx := rand.Intn(len(rngBuf))
				tet[i*6 + j].DistributeTetra(rngBuf[bufIdx], vecBuf)
				
				density := cicInterpolate(vecBuf, ghd, grid)
				
				out[i] += density / 6
			}
		}
	})
}

func ngcInterpolate(
	vecs []geom.Vec, ghd *io.GridHeader, grid []float64,
) float64 {
	L := float32(ghd.Loc.Span[0] - ghd.Loc.PixelWidth)
	dx := float32(ghd.Loc.PixelWidth)
	N := int(ghd.Loc.PixelSpan[0])
	origin := [3]float32{
		float32(ghd.Loc.Origin[0]), float32(ghd.Loc.Origin[1]),
		float32(ghd.Loc.Origin[2]),
	}

	sum := 0.0
	for _, vec := range vecs {
		for k := 0; k < 3; k++ {
			vec[k] -= origin[k]
			if vec[k] < 0 {
				vec[k] += L
			} else if vec[k] >= L {
				vec[k] -= L
			}
		}
		ix, iy, iz := int(vec[0]/dx), int(vec[1]/dx), int(vec[2]/dx)
		sum += grid[ix + iy*N + iz*N*N]
	}

	return sum / float64(len(vecs))
}

func cicInterpolate(
	vecs []geom.Vec, ghd *io.GridHeader, grid []float64,
) float64 {
	L := float32(ghd.Cosmo.BoxWidth)
	dx := float32(ghd.Loc.PixelWidth)
	N := int(ghd.Loc.PixelSpan[0])
	N2 := int(ghd.Loc.PixelSpan[0]*ghd.Loc.PixelSpan[1])
	origin := [3]float32{
		float32(ghd.Loc.Origin[0]), float32(ghd.Loc.Origin[1]),
		float32(ghd.Loc.Origin[2]),
	}

	periodic := ghd.Loc.Span[0] >= ghd.Cosmo.BoxWidth
	nTrue := int(ghd.Loc.PixelSpan[0]) - 1

	sum := 0.0
	for _, vec := range vecs {
		for k := 0; k < 3; k++ {
			vec[k] -= origin[k]
			if vec[k] < 0 {
				vec[k] += L
			} else if vec[k] >= L {
				vec[k] -= L
			}
		}
		
		xp, yp, zp := float64(vec[0]/dx), float64(vec[1]/dx), float64(vec[2]/dx)
		
		// Floor calls neccessary if xp - 0.5 wraps around.
		ix0 := int(math.Floor(xp - 0.5))
		iy0 := int(math.Floor(yp - 0.5))
		iz0 := int(math.Floor(zp - 0.5))
		
		xc, yc, zc := float64(ix0)+0.5, float64(iy0)+0.5,  float64(iz0)+0.5		
		dx, dy, dz := xp - xc, yp - yc, zp - zc
		tx, ty, tz := 1 - dx, 1 - dy, 1 - dz
		ix1, iy1, iz1 := ix0 + 1, iy0 + 1, iz0 + 1
		
		if periodic {
			if ix1 == nTrue { ix1 = 0 }
			if iy1 == nTrue { iy1 = 0 }
			if iz1 == nTrue { iz1 = 0 }
			
			if ix0 == -1 { ix0 = nTrue -1 }
			if iy0 == -1 { iy0 = nTrue -1 }
			if iz0 == -1 { iz0 = nTrue -1 }
		}
		
		sum += grid[(ix0) + (iy0)*N + (iz0)*N2]*tx*ty*tz
		sum += grid[(ix1) + (iy0)*N + (iz0)*N2]*dx*ty*tz
		sum += grid[(ix0) + (iy1)*N + (iz0)*N2]*tx*dy*tz
		sum += grid[(ix1) + (iy1)*N + (iz0)*N2]*dx*dy*tz
		sum += grid[(ix0) + (iy0)*N + (iz1)*N2]*tx*ty*dz
		sum += grid[(ix1) + (iy0)*N + (iz1)*N2]*dx*ty*dz
		sum += grid[(ix0) + (iy1)*N + (iz1)*N2]*tx*dy*dz
		sum += grid[(ix1) + (iy1)*N + (iz1)*N2]*dx*dy*dz
	}
	
	return sum / float64(len(vecs))	
}

func TetraRngBuf(pts int) []geom.Vec {
    gen := gtet_rand.NewTimeSeed(gtet_rand.Tausworthe)

	unitBuf := make([]geom.Vec, pts)
	for j := range unitBuf {
		for k := 0; k < 3; k++ {
			unitBuf[j][k] = float32(gen.Uniform(0, 1))
		}
	}
	geom.DistributeUnit(unitBuf)
	
    return unitBuf
}

// DeformationEig reutrns the eigenvalues of a deformation tensor. Requires a
// length 3 buffer of complex128 values.
func DeformationEig(def *mat.Dense, eig *mat.Eigen, buf []complex128) (
	l1, l2, l3 float64,
) {
	ok := eig.Factorize(def, mat.EigenRight)
	if !ok {
		panic(fmt.Sprintf("decomposition of %v failed", def))
	}
	val := eig.Values(buf)
	
	return sort3(real(val[0]), real(val[1]), real(val[2]))
}

func sort3(x, y, z float64) (l1, l2, l3 float64) {
	min, max := x, x
	if y > max {
		max = y
	} else if y < min {
		min = y
	}

	if z > max {
		max = z
	} else if z < min {
		min = z
	}

	return max, (x+y+z) - (min+max), min
}
