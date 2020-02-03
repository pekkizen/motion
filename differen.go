package motion
import "math"

func (c *BikeCalc) acceleration(v, power float64) (a, a1, a2 float64) {

	ov := 1.0 / v
	oMassKin := 1.0 / c.massKin
	v += c.wind
	cDrag := c.cDrag
	if v < 0 {
		cDrag = -cDrag
	}
	power *= ov
	a  = oMassKin * (power - c.fGR - cDrag*v*v) // dv/dt = acceleration
	a1 = oMassKin * (-power*ov - 2*cDrag*v)     // da/dv
	a2 = oMassKin * (power*ov*ov - cDrag) * 2   // d2a/dv2
	return 
}

// DeltaVelDist laskee (alkunopeudella v0 ja teholla power)
// ajassa Δtime tehdyn matkan Δdist ja vauhdin muutoksen Δvel.
//
func (c *BikeCalc) DeltaVelDist(Δtime, v0, power float64) (Δvel, Δdist float64) {
	const (
		n2 = 1.0 / 2
		n3 = n2  / 3
		n4 = n3  / 4
	)
	a, a1, a2 := c.acceleration(v0, power)
	// a		= 	dv/dt = acceleration
	// a1		= 	da/dv
	// a2		= 	d2a/dv2
	// jerk		=  	da/dt	=	d2v/dt2
	// jounce	=  	d2a/dt2	=	d3v/dt3
	jerk := a * a1 
	jounce := a * (a*a2 + a1*a1) 
	Δvel = Δtime*(a + Δtime*(jerk*n2 + Δtime*jounce*n3))
	Δdist = Δtime*(v0 + Δtime*(a*n2 + Δtime*(jerk*n3 + Δtime*jounce*n4)))
	return
}

// DeltaVelTime laskee (alkunopeudella v0 ja teholla power) matkaan Δdist
// kuluvan ajan Δtime ja vauhdin muutoksen Δvel.
//
func (c *BikeCalc) DeltaVelTime(Δdist, v0, power float64) (Δvel, Δtime float64) {
	const n2 = 1.0 / 2
	const n3 = n2  / 3

	a, a1, a2 := c.acceleration(v0, power)
	// a	= 	dv/dt = acceleration
	// a1	= 	da/dv
	// a2	= 	d2a/dv2

	ov := 1.0 / v0

	//f on nopeuden v 1. derivaatta matkan x suhteen: a/v (m/s/m)
	f := a * ov               //  dv/dx
	f1 := (a1 - f) * ov       //  df/dv
	f2 := (a2 - 2*f1) * ov    //  d2f/dv2
	fx1 := f * f1             //  df/dx		= d2v/dx2
	fx2 := f * (f*f2 + f1*f1) //  d2f/dx2	= d3v/dx3
	Δvel = Δdist*(f + Δdist*(fx1*n2 + Δdist*fx2*n3))

	//g on ajan t 1. derivaatta matkan x suhteen: 1/v (s/m)
	g := ov                   //  dt/dx
	g1 := -ov * ov            //  dg/dv
	g2 := -g1 * 2 * ov        //  d2g/dv2
	gx1 := f * g1             //  dg/dx		= d2t/dx2
	gx2 := f * (f*g2 + f1*g1) //  d2g/dx2	= d3t/dx3
	Δtime = Δdist*(g + Δdist*(gx1*n2 + Δdist*gx2*n3))

	//Tämä toimii melkein yhtä hyvin
	// Δtime = Δdist / (v0 + 0.5*Δvel)
	return
}

// DeltaDist laskee alkunopeudesta v0, teholla power matkan Δdist,
// joka tarvitaan nopeuden muutokseen v0 -> v1. 
//
// Ratkaistaan dist ja muut tarvittavat yhtälöistä alla. Muut tarvittavat, jotta kutsuvan
// ohjelman ei tarvitse laskea niitä toiseen kertaan,
// fRider = -(power/v0 + power/v1) / 2
// fDrag = (Fdrag(v0) + Fdrag(v1)) / 2
// cMassKin*(v0^2 - v1^2) = dist * (fGR + fRider + fDrag)
//
func (c *BikeCalc) DeltaDist(v0, v1, power float64) (Δdist, fDrag, forces float64) {
	
	Δdist = c.cMassKin * (v0*v0 - v1*v1)
	forces = c.fGR - power*(v0+v1)/(2*v0*v1)
	v0 += c.wind; v1 += c.wind
	fDrag = c.cDrag * 0.5*(v0*abs(v0) + v1*abs(v1))
	forces += fDrag 
	Δdist /= forces
	return
}

// JARRUTUSFUNKTIOT
// DistBrake(MaxEntryVel(dist, vExit), vExit) == dist
// DistBrake(VelFromBrakeDist(dist), 0) == dist
// Funktioilla VelFromBrakeDist ja MaxEntryVel lasketaan nopeusrajoituksia tyynessä.

// DistBrake laskee jarrutusmatkan nopeudesta vEntry nopeuteen vExit (tuuli mukana).
// Ratkaistaan dist ja lasketaan muut  yhtälöistä:
// fDrag = 0.5*cDrag*((vEntry+wind)^2 + (vExit+wind)^2)
// cMassKin*(vEntry^2 - vExit^2) = dist*(fGR + fBrake + fDrag)
//
func (c *BikeCalc) DistBrake(v0, v1 float64) (Δdist, fDrag, forces float64) {

	Δdist = c.cMassKin * (v0*v0 - v1*v1)
	v0 += c.wind
	v1 += c.wind
	fDrag = c.cDrag * 0.5*(v0*abs(v0) + v1*abs(v1))
	forces = c.fGR + fDrag + c.fBrake
	Δdist /= forces
	return
}

// MaxEntryVel laskee maksiminopeuden, josta etäisyydellä dist on
// mahdollista jarruttaa nopeuteen vExit (tuuli=0).
// Ratkaistaan vEntry yhtälöstä:
// cMassKin*(vEntry^2 - vExit^2) = dist*(fGR + fBrake + 0.5*cDrag*(vEntry^2 + vExit^2))
//
func (c *BikeCalc) MaxEntryVel(dist, vExit float64) float64 {

	v2 := vExit * vExit
	s := c.cMassKin - 0.5*dist*c.cDrag
	if s <= 0 {
		return maxSpeed
	}
	d := dist*(c.fGR + c.fBrake + 0.5*c.cDrag*v2) + c.cMassKin*v2
	d /= s
	if d > v2 {
		return math.Sqrt(d)
	}
	return vExit
}

// VelFromBrakeDist laskee vauhdin, josta etäisyydellä dist voidaan jarruttamalla pysähtyä (tuuli=0).
// VelFromBrakeDist(dist) == MaxEntryVel(dist, 0)
// VelFromBrakeDist on pelkästään hieman kevyempi laskea
//
// Ratkaistaan vEntry yhtälöstä:
// cMassKin*vEntry^2 = dist*(fGR + fBrake + 0.5*cDrag*vEntry^2)
//
func (c *BikeCalc) VelFromBrakeDist(dist float64) float64 {

	forces := c.fGR + c.fBrake
	if forces <= 0 {
		return 0
	}
	d := c.cMassKin - 0.5*dist*c.cDrag
	if d <= 0 {
		return maxSpeed
	}
	return math.Sqrt(dist * forces / d)
}
