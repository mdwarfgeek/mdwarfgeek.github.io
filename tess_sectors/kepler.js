const TWOPI = 2.0*Math.PI;

const KEPLER_PREC    = 1.0e-13;
const KEPLER_MAXITER = 100;

function kepler (ma, ecc) {
  var ea;

  // Reduce mean anomaly to [-pi, pi].  Peculiar semantics of Javascript's
  // round function don't really matter here.
  var tmp = ma / TWOPI;
  ma = TWOPI * (tmp - Math.round(tmp));

  // Eccentric or circular?
  if(ecc > 0) {
    // For eccentric orbits, use cubic approximation from Mikkola (1987)
    // to provide initial guess of eccentric anomaly.
    
    // Eq. 9a
    var tmp = 1.0 / (4 * ecc + 0.5);
    
    var alpha = (1.0 - ecc) * tmp;
    var beta = 0.5 * ma * tmp;
    
    // Eq. 9b
    tmp = Math.sqrt(beta*beta + alpha*alpha*alpha);
    if(beta < 0)
      tmp = -tmp;

    var z = Math.cbrt(beta + tmp);
    
    // Eq. 9c: initial value of sin(E/3)
    var ste = z - alpha / z;
    
    // Eq. 7: 5th order correction term
    ste -= 0.078 * ste*ste*ste*ste*ste / (1.0 + ecc);
    
    // Eq. 8: eccentric anomaly
    ea = ma + ecc * ste * (3.0 - 4.0 * ste*ste);
    
    // Refine solution of Kepler's equation using Newton's method
    for(i = 0; i < KEPLER_MAXITER; i++) {
      se = Math.sin(ea);
      ce = Math.cos(ea);

      var f = ea - ecc * se - ma;
      var df = 1.0 - ecc * ce;
      var delta = f / df;
      
      ea -= delta;
      
      if(Math.abs(delta) < KEPLER_PREC) {
	// I think that's enough...
        break;
      }
    }
    
    if(i >= KEPLER_MAXITER)
      console.log("kepler: iteration limit reached");
  }
  else
    ea = ma;

  return ea;
}
