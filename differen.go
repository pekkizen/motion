package motion

import "math"

func (c *BikeCalc) Acceleration(v, power float64) float64 {
	// 1. derivative of speed v with respect to time t: dv/dt = acceleration a = F/m.
	// where m is mass and F is sum of resisting/assisting forces.

	power /= v
	v += c.wind
	return (power - c.fGR - c.cDrag*v*math.Abs(v)) * c.oMassKin
}

// At the end of this file are these two functions, DeltaTime and DeltaDist
// implemented by 3 derivatives of speed, respect to time and distance.
// Three derivatives are not enough for bigger delta time and distance and
// using more derivatives is not practical. Also single derivative is not enough,
// although mostly not worse than 3 derivatives, but estimating "average/middle"
// acceleration below seems to work very well.

func (c *BikeCalc) DeltaTime(Δtime, v, power float64) (Δvel, Δdist float64) {

	Δvel = Δtime * c.Acceleration(v, power)
	Δvel = Δtime * c.Acceleration(v+0.5*Δvel, power)
	Δdist = Δtime * (v + 0.5*Δvel)
	return
}

func (c *BikeCalc) DeltaDist(Δdist, v, power float64) (Δvel, Δtime float64) {
	// 1. derivative of speed v with respect to distance s:  dv/ds = a/v,
	// where a is acceleration.

	Δvel = Δdist * c.Acceleration(v, power) / v
	w := v + 0.5*Δvel
	Δvel = Δdist * c.Acceleration(w, power) / w
	Δtime = Δdist / (v + 0.5*Δvel)
	return
}

// DeltaVel laskee alkunopeudesta v0, teholla power matkan Δdist,
// joka tarvitaan nopeuden muutokseen v0 -> v0 + Δvel.
//
// Change in kinetic energy ΔKE = distance x force
// fGR 		= gravity and rolling resistance force
// fRider 	= -(power/v0 + power/v1) / 2
// fDrag 	= (Fdrag(v0) + Fdrag(v1)) / 2
// ΔKE 		= 0.5 * mass * (v0^2 - v1^2)
// ΔKE		= Δdist * (fGR + fRider + fDrag)
// ΔKE/mass	= Δdist * acceleration
func (c *BikeCalc) DeltaVel(Δvel, v0, power float64) (Δdist, Δtime, fDrag, fSum float64) {

	v1 := v0 + Δvel
	ΔKE := 0.5 * c.massKin * (v0*v0 - v1*v1)
	fRider := -power * (v0 + v1) / (2 * v0 * v1)
	v0 += c.wind
	v1 += c.wind
	fDrag = c.cDrag * 0.5 * (v0*math.Abs(v0) + v1*math.Abs(v1))
	fSum = c.fGR + fRider + fDrag
	Δdist = ΔKE / fSum // forces near 0 must be checked in calling program
	Δtime = -Δvel * c.massKin / fSum
	// Δtime = Δdist / (v0 + 0.5*Δvel) 		// same to 7 digits ride time
	return
}

// JARRUTUSFUNKTIOT
// DistBrake(MaxEntryVel(dist, vExit), vExit) == dist
// DistBrake(VelFromBrakeDist(dist), 0) == dist
// Funktioilla VelFromBrakeDist ja MaxEntryVel lasketaan nopeusrajoituksia tyynessä.

// DistBrake laskee jarrutusmatkan nopeudesta vEntry nopeuteen vExit (tuuli mukana).
// Ratkaistaan dist ja lasketaan muut  yhtälöistä:
// fDrag = 0.5*cDrag*((vEntry+wind)^2 + (vExit+wind)^2)
// 0.5 * mass * (vEntry^2 - vExit^2) = dist*(fGR + fBrake + fDrag)
func (c *BikeCalc) DistBrake(v0, v1 float64) (Δdist, fDrag, fSum float64) {

	ΔKE := 0.5 * c.massKin * (v0*v0 - v1*v1)
	v0 += c.wind
	v1 += c.wind
	fDrag = c.cDrag * 0.5 * (v0*math.Abs(v0) + v1*math.Abs(v1))
	fSum = c.fGR + fDrag + c.fBrake
	Δdist = ΔKE / fSum
	return
}

// MaxEntryVel laskee maksiminopeuden, josta etäisyydellä dist on
// mahdollista jarruttaa nopeuteen vExit (tuuli=0).
// Solve vEntry from equation:
// 0.5 * mass  * (vEntry^2 - vExit^2) = dist*(fGR + fBrake + 0.5*cDrag*(vEntry^2 + vExit^2))
func (c *BikeCalc) MaxEntryVel(dist, vExit float64) float64 {

	w := vExit * vExit
	s := 0.5 * (c.massKin - dist*c.cDrag)
	if s <= 0 {
		return maxSpeed
	}
	d := dist*(c.fGR+c.fBrake+0.5*c.cDrag*w) + 0.5*c.massKin*w
	d /= s
	if d > w {
		return math.Sqrt(d)
	}
	return vExit
}

// VelFromBrakeDist laskee vauhdin, josta etäisyydellä dist voidaan jarruttamalla pysähtyä (tuuli=0).
// VelFromBrakeDist(dist) == MaxEntryVel(dist, 0)
//
// Solve vEntry from equation::
// 0.5 * mass  * vEntry^2 = dist*(fGR + fBrake + 0.5*cDrag*vEntry^2)
func (c *BikeCalc) VelFromBrakeDist(dist float64) float64 {

	fSum := c.fGR + c.fBrake
	if fSum <= 0 {
		return 0
	}
	d := 0.5 * (c.massKin - dist*c.cDrag)
	if d <= 0 {
		return maxSpeed
	}
	return math.Sqrt(dist * fSum / d)
}

// These below are not used any more
func (c *BikeCalc) accelerationDerivatives(v, power float64) (a, a1, a2 float64) {

	ov := 1.0 / v
	omk := 1.0 / c.massKin
	v += c.wind
	cDrag := c.cDrag
	if v < 0 {
		cDrag = -cDrag
	}
	power *= ov
	a = omk * (power - c.fGR - cDrag*v*v) // dv/dt = acceleration a = F/m
	a1 = omk * (-power*ov - 2*cDrag*v)    // da/dv
	a2 = omk * (power*ov*ov - cDrag) * 2  // d2a/dv2
	return
}

func (c *BikeCalc) DeltaDist3(Δdist, v0, power float64) (Δvel, Δtime float64) {
	const n2 = 1.0 / 2
	const n3 = n2 / 3
	ov := 1.0 / v0

	a, a1, a2 := c.accelerationDerivatives(v0, power)
	// a	= 	dv/dt = acceleration
	// a1	= 	da/dv
	// a2	= 	d2a/dv2
	// f1 on nopeuden v0 1. derivaatta matkan x suhteen: a/v (m/s/m)

	f1 := a * ov               //  dv/dx
	b1 := (a1 - f1) * ov       //  df/dv
	b2 := (a2 - 2*b1) * ov     //  d2f/dv2
	f2 := f1 * b1              //  df/dx	= d2v/dx2
	f3 := f1 * (f1*b2 + b1*b1) //  d2f/dx2	= d3v/dx3
	Δvel = Δdist * (f1 + Δdist*(f2*n2+Δdist*f3*n3))
	Δtime = Δdist / (v0 + 0.5*Δvel)

	// g on ajan t 1. derivaatta matkan x suhteen: 1/v (s/m)
	// g := ov                   	//  dt/dx
	// g1 := -ov * ov            	//  dg/dv
	// g2 := -g1 * 2 * ov        	//  d2g/dv2
	// gx1 := f1 * g1             	//  dg/dx	= d2t/dx2
	// gx2 := f1 * (f1*g2 + b1*g1) 	//  d2g/dx2	= d3t/dx3
	// Δtime = Δdist * (g + Δdist*(gx1*n2+Δdist*gx2*n3))
	return
}

func (c *BikeCalc) DeltaTime3(Δtime, v0, power float64) (Δvel, Δdist float64) {
	const (
		n2 = 1.0 / 2
		n3 = n2 / 3
		n4 = n3 / 4
	)
	a1, b1, b2 := c.accelerationDerivatives(v0, power)
	// a1		= 	dv/dt = acceleration
	// b1		= 	da/dv
	// b2		= 	d2a/dv2
	// a2		=  	da/dt	=	d2v/dt2
	// a3 		=  	d2a/dt2	=	d3v/dt3
	a2 := a1 * b1
	a3 := a1 * (a1*b2 + b1*b1)
	Δvel = Δtime * (a1 + Δtime*(a2*n2+Δtime*a3*n3))
	Δdist = Δtime * (v0 + 0.5*Δvel)
	// Δdist = Δtime * (v0 + Δtime*(a1*n2+Δtime*(a2*n3+Δtime*a3*n4)))
	return
}
