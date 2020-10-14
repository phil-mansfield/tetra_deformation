package lib

import (
	"gonum.org/v1/gonum/mat"
	"github.com/phil-mansfield/gotetra/render/geom"
	"github.com/phil-mansfield/gotetra/render/io"
)

const (
	AllTetra int = iota
)

// Order specifies whcih tetrahedra will actually be analyzed and the order that
// they will be analyzed in. Currently, this only means specifying a skip.
type Order struct {
	Skip int // Must be a power of 2.
	TetraType int // Int indicating which tetrahedra to create for each point.
	// (This is currently unused)
}

// Buffers returns the vector buffer and tetra buffer that would be needed to
// read tetrahedra from fname in a given order.
func Buffers(order *Order, fname string) (
	[]geom.Vec, []geom.Tetra, []*mat.Dense,
) {
	hd := &io.SheetHeader{ }
	err := io.ReadSheetHeaderAt(fname, hd)
	if err != nil { panic(err.Error()) }

	x := make([]geom.Vec, hd.GridCount)

	tetWidth := int(hd.SegmentWidth) / order.Skip
	tet := make([]geom.Tetra, 6 * tetWidth*tetWidth*tetWidth)

	def := make([]*mat.Dense, tetWidth*tetWidth*tetWidth)
	for i := range def {
		def[i] = mat.NewDense(3, 3, make([]float64, 9))
	}
	
	return x, tet, def
}

func Read(
	order *Order, fname string,
	x []geom.Vec, tet []geom.Tetra, def []*mat.Dense,
) {
	hd := &io.SheetHeader{ }
	err := io.ReadSheetHeaderAt(fname, hd)
	if err != nil { panic(err.Error()) }

	err = io.ReadSheetPositionsAt(fname, x)
	if err != nil { panic(err.Error()) }

	Ns, Ng, skip := int(hd.SegmentWidth), int(hd.GridWidth), int(order.Skip)
	dX := float64(order.Skip) * hd.TotalWidth / float64(hd.CountWidth) 
	
	j := 0
	idxBuf := &geom.TetraIdxs{ }
	vecBuf := make([]geom.Vec, 4)
	defIdx := make([]int, 4)
	for iz := 0; iz < Ns; iz += skip {
		for iy := 0; iy < Ns; iy += skip {
			for ix := 0; ix < Ns; ix += skip {
				i := ix + iy*Ng + iz*Ng*Ng

				initDefIdx(i, Ng, order.Skip, defIdx)
				initDef(def[j/6], defIdx, x, dX, hd.TotalWidth)
				
				for dir := 0; dir < 6; dir++ {
					idxBuf.Init(int64(i), int64(Ng), int64(order.Skip), dir)
					for i := range idxBuf {
						vecBuf[i] = x[idxBuf[i]]
					}
					
					boundVecBuf(vecBuf, float32(hd.TotalWidth))
					tet[j].Init(&vecBuf[0], &vecBuf[1], &vecBuf[2], &vecBuf[3])

					j++
				}
			}
		}
	}
}

// boundVecBuf
func boundVecBuf(buf []geom.Vec, L float32) {
	L2 := L / 2
	
	for k := 0; k < 3; k++ {
		for i := 1; i < len(buf); i++ {
			if dx := buf[i][k] - buf[0][k]; dx > L2 {
				buf[i][k] -= L
			} else if dx < -L/2 {
				buf[i][k] += L
			}
		}
	}
}

func initDefIdx(i, Ng, skip int, buf []int) {
	buf[0] = i
	buf[1] = i + Ng*Ng*skip // x
	buf[2] = i + Ng*skip // y
	buf[3] = i + skip // x
}

func initDef(def *mat.Dense, idx []int, x []geom.Vec, dX, L float64) {
	for i := 0; i < 3; i++ {
		for j := 0; j < 3; j++ {
			dx := float64(x[idx[i + 1]][j] - x[idx[0]][j])
			if i == j { dx -= dX }
			if dx > L/2 {
				dx -= L
			} else if dx < -L/2 {
				dx += L
			}
			def.Set(i, j, dx / dX)
		}
	}
}
