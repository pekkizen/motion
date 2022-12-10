package motion

import "math"

func (c *BikeCalc) Acceleration(v, power float64) float64 {
	// 1. derivative of speed v with respect to time t: dv/dt = acceleration a = F/m.
	// where m is mass and F is sum of resisting/assisting forces.
	power /= v
	v += c.wind
	return (power - c.fGR - c.cDrag*math.Abs(v)*v) * c.oMassKin
}

func (c *BikeCalc) accelerationDerivatives(v, power float64) (dv_dt, da_dv, d2a_dv2 float64) {
	ov := 1.0 / v
	v += c.wind
	cDrag := math.Copysign(c.cDrag, v)
	power *= ov
	dv_dt = c.oMassKin * (power - c.fGR - cDrag*v*v) // dv/dt = acceleration a = F/m
	da_dv = c.oMassKin * (-power*ov - 2*cDrag*v)     // da/dv
	d2a_dv2 = c.oMassKin * (power*ov*ov - cDrag) * 2 // d2a/dv2
	return
}

func (c *BikeCalc) DeltaTime(Δtime, v0, power float64) (Δvel, Δdist float64) {

	a, b1, b2 := c.accelerationDerivatives(v0, power)
	// a1		= 	dv/dt = acceleration
	// b1		= 	da/dv
	// b2		= 	d2a/dv2
	// a2		=  	da/dt	=	d2v/dt2
	// a3 		=  	d2a/dt2	=	d3v/dt3
	a2 := a * b1
	a3 := a * (a*b2 + b1*b1)
	Δvel = Δtime * (a + Δtime*(0.5*a2+(1.0/6)*Δtime*a3))
	Δdist = Δtime * (v0 + 0.5*Δvel)
	return
}

func (c *BikeCalc) DeltaDist(Δdist, v0, power float64) (Δvel, Δtime float64) {
	ov := 1.0 / v0
	a, a1, a2 := c.accelerationDerivatives(v0, power)

	// f is the 1. derivative of speed v with respect to distance x:  dv/dx = a/v
	f := a * ov              //  dv/dx
	b1 := (a1 - f) * ov      //  df/dv
	b2 := (a2 - 2*b1) * ov   //  d2f/dv2
	f2 := f * b1             //  df/dx	= d2v/dx2
	f3 := f * (f*b2 + b1*b1) //  d2f/dx2	= d3v/dx3

	Δvel = Δdist * (f + Δdist*(0.5*f2+(1.0/6)*Δdist*f3))
	Δtime = Δdist / (v0 + 0.5*Δvel)
	return
}

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
	Δdist = ΔKE / fSum 						// check fSum for 0 in calling program
	Δtime = -Δvel * c.massKin / fSum		// fSUM = 0 => acceleration = 0
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
