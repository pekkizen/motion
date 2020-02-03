package motion

import (
	"testing"
	"math"
	// "math/rand"
	rand "prng"
)

var SinkTest float64

var Cal = Calculator()



func (c *BikeCalc) functionTst(v float64) float64 {
	c.callsFunc++
	return v*(c.fGR + c.cDrag * abs(v+c.wind)*(v+c.wind)) - c.power
}
func (c *BikeCalc) functionTst2(v float64) float64 {
	c.callsFunc++
	return v*(c.fGR + c.cDrag * abs2(v+c.wind)*(v+c.wind)) - c.power
	// return v*(c.fGR + c.cDrag * signedSquare(v+c.wind)) - c.power
}

func signedSquare(x float64) float64 {
	return math.Float64frombits(math.Float64bits(x*x) | (math.Float64bits(x) & (1 << 63)))
}
var rng = rand.New(0)
// var rng = rand.New(rand.NewSource(99))

func BenchmarkRand(b *testing.B) {
	sign := 1.0
	y :=  0.0
	for n := 0; n < b.N; n++ {
		x := rand.Float64() * sign
		sign *= -1
		y += x
	}
	SinkTest = y
}
func BenchmarkNull(b *testing.B) {	
	var k int
	for n := 0; n < b.N; n++ {
		k += n
	}
	SinkTest += float64(k)
}

func BenchmarkAbs(b *testing.B) {
	sign := 100.0
	y := 0.0
	for n := 0; n < b.N; n++ {
		x := rand.Float64() * sign
		sign *= -1
		// y += abs(x)
		y += abs2(x)
		y += abs2(x)
		y += abs2(x)
		y += abs2(x)
		y += abs2(x)
		y += abs2(x)
		y += abs2(x)
		y += abs2(x)
		y += abs2(x)
		y += abs2(x)
	}
	SinkTest = y
}

func BenchmarkSignedSquare(b *testing.B) {
	sign := 100.0
	y := 10.0
	for n := 0; n < b.N; n++ {
		x := rand.Float64() * sign
		sign *= -1
		y += signedSquare(x)
	}
	SinkTest += y
}
func BenchmarkAbsXx(b *testing.B) {
	sign := 100.0
	y := 0.0
	for n := 0; n < b.N; n++ {
		x := rand.Float64() * sign
		// sign *= -1
		y += x * abs2(x)
		// y += x * abs(x)
	}
	SinkTest += y
}


func BenchmarkFunction(b *testing.B) {
	Cal.SetWeight(100)
	Cal.SetGrade(0.05)
	Cal.SetWind(-3)
	Cal.SetPower(100)
	y := 0.0
	for n := 0; n < b.N; n++ {
		y += Cal.functionTst(5) 
	}
	SinkTest = y
}

func BenchmarkSetGrade(b *testing.B) {
	Cal.SetWeight(100)
	for n := 0; n < b.N; n++ {
		Cal.SetGrade(0.01)
	}
}
func BenchmarkSetGradeExact(b *testing.B) {
	Cal.SetWeight(100)
	for n := 0; n < b.N; n++ {
		Cal.SetGradeExact(0.01)
	}
}
func BenchmarkNewtonRaphson(b *testing.B) {
	// var newton =  (*BikeCalc).NewtonRaphson
	Cal.SetWeight(100)
	Cal.SetGrade(0.0001)
	Cal.SetWind(0)
	Cal.SetCrr(0.007)
	Cal.SetCdA(0.7)
	v := 5.0
	for n := 0; n < b.N; n++ {
		v, _ = Cal.NewtonRaphson(100, 0.1, v)
		v += 0.25
	}
	SinkTest = v
}
func BenchmarkNewtonRaphsonHalley(b *testing.B) {
	Cal.SetWeight(100)
	Cal.SetGrade(0.0001)
	Cal.SetWind(-1)
	Cal.SetCrr(0.007)
	Cal.SetCdA(0.7)
	v := 5.0
	for n := 0; n < b.N; n++ {
		// v = float64(n % 9 + 1)
		v, _ = Cal.NewtonRaphsonHalley(100, 0.1, v)
		if n % 2 == 0 {
			v += 0.35
		}
		
	}
	SinkTest = v
}
func BenchmarkSingleQuadratic(b *testing.B) {
	Cal.SetWeight(100)
	Cal.SetGrade(0.0)
	Cal.SetWind(0)
	Cal.SetCrr(0.007)
	Cal.SetCdA(0.7)
	for n := 0; n < b.N; n++ {
		Cal.SingleQuadratic(100, 1.5, 5)
	}
}
func BenchmarkSingleQuadraticTriple(b *testing.B) {
	Cal.SetWeight(100)
	Cal.SetGrade(0.0)
	Cal.SetWind(0)
	Cal.SetCrr(0.007)
	Cal.SetCdA(0.7)
	for n := 0; n < b.N; n++ {
		Cal.SingleQuadraticTriple(100, 1.5, 5)
	}
}
func BenchmarkBracket(b *testing.B) {
	Cal.SetWeight(100)
	Cal.SetGrade(0.0)
	Cal.SetPower(100)
	Cal.SetWind(0)
	Cal.SetCrr(0.007)
	Cal.SetCdA(0.7)
	for n := 0; n < b.N; n++ {
		Cal.Bracket(2.5, 5)
	}
}
func BenchmarkDealtaDist(b *testing.B) {
	Cal.SetWeight(100)
	Cal.SetGrade(0.0)
	Cal.SetWind(0)
	Cal.SetCrr(0.007)
	Cal.SetCdA(0.7)
	dist := 0.0
	for n := 0; n < b.N; n++ {
		dist, _, _ = Cal.DeltaDist(3, 6, 200)
	}
	SinkTest = dist
}