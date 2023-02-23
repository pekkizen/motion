package motion

import (
	"math"
	"testing"

	// "math/rand"
	rand "prng"
)

var fsink float64

var c = Calculator()

// func (c *BikeCalc) functionTst(v float64) float64 {
// 	c.callsFunc++
// 	return v*(c.fGR + c.cDrag * math.Abs(v+c.wind)*(v+c.wind)) - c.power
// }

// func signedSquare(x float64) float64 {
// 	return math.Float64frombits(math.Float64bits(x*x) | (math.Float64bits(x) & (1 << 63)))
// }

// var rng = rand.New(0)
// var rng = rand.New(rand.NewSource(99))

func fastExp(x float64) float64 {
	const o4096 = 1.0 / 4096.0
	x = 1 + x*o4096
	x *= x
	x *= x
	x *= x
	x *= x
	x *= x
	x *= x
	x *= x
	x *= x
	x *= x
	x *= x
	x *= x
	x *= x
	return x
}
func BenchmarkFastExp(b *testing.B) {

	y := 10.0
	for n := 0; n < b.N; n++ {
		y += fastExp(1)
	}
	fsink += y
}
func BenchmarkExp(b *testing.B) {
	y := 10.0
	for n := 0; n < b.N; n++ {
		y += math.Exp(1)
	}
	fsink += y
}

func BenchmarkRand(b *testing.B) {
	sign := 1.0
	y := 0.0
	for n := 0; n < b.N; n++ {
		x := rand.Float64() * sign
		sign *= -1
		y += x
	}
	fsink = y
}
func BenchmarkNull(b *testing.B) {
	var k int
	for n := 0; n < b.N; n++ {
		k += n
	}
	fsink += float64(k)
}

func BenchmarkSignedSquare(b *testing.B) {
	y := 10.0
	z := 1.0
	x := 10.0
	for n := 0; n < b.N; n++ {
		x++
		// x *= -1
		y = signedSquare(x+z)
		// y = (x + z) * math.Abs(x+z)
	}
	fsink = y
}

func BenchmarkFunction(b *testing.B) {
	c.SetWeight(100)
	c.SetGrade(0.05)
	c.SetWind(1)
	c.SetPower(100)
	v := 10.0
	y := 0.0
	for n := 0; n < b.N; n++ {
		// v = float64(n)
		// v++
		y = c.function(v)
		// y = v*(c.fGR+c.cDrag*math.Abs(v+c.wind)*(v+c.wind)) - c.power
	
	}
	fsink = y
}

func BenchmarkSetGrade(b *testing.B) {
	c.SetWeight(100)
	for n := 0; n < b.N; n++ {
		c.SetGrade(0.01)
	}
}

func BenchmarkNewtonRaphson(b *testing.B) {
	// var newton =  (*BikeCalc).NewtonRaphson
	c.SetWeight(100)
	c.SetGrade(-0.02)
	c.SetWind(0)
	c.SetCrr(0.007)
	c.SetCdA(0.7)
	v := 5.0
	p := 150.0
	for n := 0; n < b.N; n++ {
		v, _ = c.NewtonRaphson(p, 0.05, v)
		v += 1
	}
	fsink = v
}

func BenchmarkSingleQuadratic(b *testing.B) {
	c.SetWeight(100)
	c.SetGrade(-0.02)
	c.SetWind(0)
	c.SetCrr(0.007)
	c.SetCdA(0.7)
	v := 5.0
	p := 150.0
	for n := 0; n < b.N; n++ {
		v, _ = c.Quadratic(p, 0.5, v)
		v += 1
	}
	fsink = v
}

func BenchmarkLinear(b *testing.B) {
	c.SetWeight(100)
	c.SetGrade(-0.02)
	c.SetWind(0)
	c.SetCrr(0.007)
	c.SetCdA(0.7)
	v := 5.0
	p := 150.0
	for n := 0; n < b.N; n++ {
			v, _ = c.SingleLinear(p, 0.5, v)
		v += 1
	}
	fsink = v
}

func BenchmarkDoubleLinear(b *testing.B) {
	c.SetWeight(100)
	c.SetGrade(-0.02)
	c.SetWind(0)
	c.SetCrr(0.007)
	c.SetCdA(0.7)
	v := 5.0
	for n := 0; n < b.N; n++ {
		p := 150.0
		v, _ = c.DoubleLinear(p, 0.5, v)
		v += 1
	}
	fsink = v
}

func BenchmarkBracket(b *testing.B) {
	c.SetWeight(100)
	c.SetGrade(-0.02)
	c.SetPower(100)
	c.SetWind(0)
	c.SetCrr(0.007)
	c.SetCdA(0.7)
	v := 5.0
	for n := 0; n < b.N; n++ {
		v, _, _, _ = c.bracket(2.5, v)
		v += 5
	}
	fsink = v
}

func BenchmarkDeltaVel(b *testing.B) {
	c.SetWeight(100)
	c.SetGrade(0.0)
	c.SetWind(0)
	c.SetCrr(0.007)
	c.SetCdA(0.7)
	dist := 0.0
	time := 0.0
	// z := 0.0
	for n := 0; n < b.N; n++ {
		// dist, time, _,_ = c.DeltaVelBrake(1, 6)
		// dist,time, _,_ = c.DeltaVelBrakeMD(1, 6)
		// dist,time, _,_ = c.DeltaVel(1, 6, 100)
		// dist,time, _,_ = c.DeltaVelMD(1, 6, 100)
		dist, time = c.DeltaTime(1, 6, 100)
		// dist, time = c.DeltaDist(1, 6, 100)
		// z = dist
	}
	fsink = dist + time
}

func BenchmarkMeanFdrag(b *testing.B) {
	c.SetWeight(100)
	c.SetGrade(0.0)
	c.SetWind(0)
	c.SetCrr(0.007)
	c.SetCdA(0.7)
	v0 := 5.0
	dVel := 2.0
	z := 0.0
	for n := 0; n < b.N; n++ {
		z = c.MeanFdrag(dVel, v0)
		// z = MeanFrider(dVel, v0 , 100)
		// z = c.cDrag * signedSquare(v0+ 0.5*dVel + c.wind)
		// z = (1.0 / 3) *  c.cDrag * (3 * v0 * (v0 +dVel) + dVel * dVel)
	}
	fsink = z 
}

