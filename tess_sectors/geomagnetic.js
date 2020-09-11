// Geomagnetic field, WMM2020.
const wmm_coef = [
//       g        h     g'     h'      n  m
  -29404.5,     0.0,   6.7,   0.0, //  1  0
   -1450.7,  4652.9,   7.7, -25.1, //  1  1
   -2500.0,     0.0, -11.5,   0.0, //  2  0
    2982.0, -2991.6,  -7.1, -30.2, //  2  1
    1676.8,  -734.8,  -2.2, -23.9, //  2  2
    1363.9,     0.0,   2.8,   0.0, //  3  0
   -2381.0,   -82.2,  -6.2,   5.7, //  3  1
    1236.2,   241.8,   3.4,  -1.0, //  3  2
     525.7,  -542.9, -12.2,   1.1, //  3  3
     903.1,     0.0,  -1.1,   0.0, //  4  0
     809.4,   282.0,  -1.6,   0.2, //  4  1
      86.2,  -158.4,  -6.0,   6.9, //  4  2
    -309.4,   199.8,   5.4,   3.7, //  4  3
      47.9,  -350.1,  -5.5,  -5.6, //  4  4
    -234.4,     0.0,  -0.3,   0.0, //  5  0
     363.1,    47.7,   0.6,   0.1, //  5  1
     187.8,   208.4,  -0.7,   2.5, //  5  2
    -140.7,  -121.3,   0.1,  -0.9, //  5  3
    -151.2,    32.2,   1.2,   3.0, //  5  4
      13.7,    99.1,   1.0,   0.5, //  5  5
      65.9,     0.0,  -0.6,   0.0, //  6  0
      65.6,   -19.1,  -0.4,   0.1, //  6  1
      73.0,    25.0,   0.5,  -1.8, //  6  2
    -121.5,    52.7,   1.4,  -1.4, //  6  3
     -36.2,   -64.4,  -1.4,   0.9, //  6  4
      13.5,     9.0,   0.0,   0.1, //  6  5
     -64.7,    68.1,   0.8,   1.0, //  6  6
      80.6,     0.0,  -0.1,   0.0, //  7  0
     -76.8,   -51.4,  -0.3,   0.5, //  7  1
      -8.3,   -16.8,  -0.1,   0.6, //  7  2
      56.5,     2.3,   0.7,  -0.7, //  7  3
      15.8,    23.5,   0.2,  -0.2, //  7  4
       6.4,    -2.2,  -0.5,  -1.2, //  7  5
      -7.2,   -27.2,  -0.8,   0.2, //  7  6
       9.8,    -1.9,   1.0,   0.3, //  7  7
      23.6,     0.0,  -0.1,   0.0, //  8  0
       9.8,     8.4,   0.1,  -0.3, //  8  1
     -17.5,   -15.3,  -0.1,   0.7, //  8  2
      -0.4,    12.8,   0.5,  -0.2, //  8  3
     -21.1,   -11.8,  -0.1,   0.5, //  8  4
      15.3,    14.9,   0.4,  -0.3, //  8  5
      13.7,     3.6,   0.5,  -0.5, //  8  6
     -16.5,    -6.9,   0.0,   0.4, //  8  7
      -0.3,     2.8,   0.4,   0.1, //  8  8
       5.0,     0.0,  -0.1,   0.0, //  9  0
       8.2,   -23.3,  -0.2,  -0.3, //  9  1
       2.9,    11.1,   0.0,   0.2, //  9  2
      -1.4,     9.8,   0.4,  -0.4, //  9  3
      -1.1,    -5.1,  -0.3,   0.4, //  9  4
     -13.3,    -6.2,   0.0,   0.1, //  9  5
       1.1,     7.8,   0.3,   0.0, //  9  6
       8.9,     0.4,   0.0,  -0.2, //  9  7
      -9.3,    -1.5,   0.0,   0.5, //  9  8
     -11.9,     9.7,  -0.4,   0.2, //  9  9
      -1.9,     0.0,   0.0,   0.0, // 10  0
      -6.2,     3.4,   0.0,   0.0, // 10  1
      -0.1,    -0.2,   0.0,   0.1, // 10  2
       1.7,     3.5,   0.2,  -0.3, // 10  3
      -0.9,     4.8,  -0.1,   0.1, // 10  4
       0.6,    -8.6,  -0.2,  -0.2, // 10  5
      -0.9,    -0.1,   0.0,   0.1, // 10  6
       1.9,    -4.2,  -0.1,   0.0, // 10  7
       1.4,    -3.4,  -0.2,  -0.1, // 10  8
      -2.4,    -0.1,  -0.1,   0.2, // 10  9
      -3.9,    -8.8,   0.0,   0.0, // 10 10
       3.0,     0.0,   0.0,   0.0, // 11  0
      -1.4,     0.0,  -0.1,   0.0, // 11  1
      -2.5,     2.6,   0.0,   0.1, // 11  2
       2.4,    -0.5,   0.0,   0.0, // 11  3
      -0.9,    -0.4,   0.0,   0.2, // 11  4
       0.3,     0.6,  -0.1,   0.0, // 11  5
      -0.7,    -0.2,   0.0,   0.0, // 11  6
      -0.1,    -1.7,   0.0,   0.1, // 11  7
       1.4,    -1.6,  -0.1,   0.0, // 11  8
      -0.6,    -3.0,  -0.1,  -0.1, // 11  9
       0.2,    -2.0,  -0.1,   0.0, // 11 10
       3.1,    -2.6,  -0.1,   0.0, // 11 11
      -2.0,     0.0,   0.0,   0.0, // 12  0
      -0.1,    -1.2,   0.0,   0.0, // 12  1
       0.5,     0.5,   0.0,   0.0, // 12  2
       1.3,     1.3,   0.0,  -0.1, // 12  3
      -1.2,    -1.8,   0.0,   0.1, // 12  4
       0.7,     0.1,   0.0,   0.0, // 12  5
       0.3,     0.7,   0.0,   0.0, // 12  6
       0.5,    -0.1,   0.0,   0.0, // 12  7
      -0.2,     0.6,   0.0,   0.1, // 12  8
      -0.5,     0.2,   0.0,   0.0, // 12  9
       0.1,    -0.9,   0.0,   0.0, // 12 10
      -1.1,     0.0,   0.0,   0.0, // 12 11
      -0.3,     0.5,  -0.1,  -0.1  // 12 12
];

const wmm_nmax = 12;
const wmm_coefmax = (wmm_nmax+1)*(wmm_nmax+2)/2;

// WGS84 constants.
const AEARTH = 6378137.0;  // m
const FEARTH = 1.0/298.257223563;

function geomagnetic (mjd, longitude, latitude, height) {
  // Transform epoch.
  var epoch = 2000.0 + (mjd - 51544.5) / 365.25;
  var t = (epoch - 2020.0);

  // sin, cos of longitude and latitude.
  var sinlam = Math.sin(longitude);
  var coslam = Math.cos(longitude);
  
  var sinphi = Math.sin(latitude);
  var cosphi = Math.cos(latitude);
  
  // Square of compression factor: b^2/a^2 = (1 - f)^2
  var cfsq = 1.0 - FEARTH;
  cfsq *= cfsq;

  // Distance from the centre * b/a.
  var n = AEARTH / Math.sqrt(cosphi*cosphi + cfsq * sinphi*sinphi);

  // Distance from Earth axis, m */
  var gu = (n      + height) * cosphi;
  var gz = (n*cfsq + height) * sinphi;

  // Normalized.
  var r = Math.hypot(gu, gz);

  var cospp = gu / r;  // cos of colatitude
  var sinpp = gz / r;  // sin of colatitude

  // Compute array of sin,cos of n*longitude.
  var sl = new Array(wmm_nmax+1);
  var cl = new Array(wmm_nmax+1);

  sl[0] = 0;
  cl[0] = 1;
  
  for(var n = 1; n <= wmm_nmax; n++) {
    sl[n] = sl[n-1] * coslam + cl[n-1] * sinlam;
    cl[n] = cl[n-1] * coslam - sl[n-1] * sinlam;
  }

  // Calculate Schmidt semi-normalized associated Legendre function p
  // and its derivative q for each (n,m) using recurrence relations
  // for argument x = sinpp.  These use n-1 and n-2 results so we need
  // to calculate n=0 and 1 explicitly before starting the loop.
  // See WMM2020_Report.pdf Eq. 4,5 for definition of these functions.
  var p = new Array(wmm_coefmax);
  var q = new Array(wmm_coefmax);

  p[0] =    1.0; // (n,m) = 0,0
  q[0] =    0.0;

  p[1] =  sinpp; // (n,m) = 1,0
  q[1] =  cospp;
  
  p[2] =  cospp; // (n,m) = 1,1
  q[2] = -sinpp; // the factor sqrt(2*(n-m)/(n+m)) = 1 here
  
  // Loop for n=2 onwards.
  var ic = 3;
  
  for(var n = 2; n <= wmm_nmax; n++) {
    var fn = n;
    var gn = n-1;
    
    for(var m = 0; m < n; m++) {
      // sqrt(n^2-m^2) p(n,m) = ((2n-1) x p(n-1,m) - sqrt((n-1)^2-m^2) p(n-2,m))
      var gmm = m*m;
      var denom = Math.sqrt(fn*fn - gmm);
      var ca = (fn + gn) / denom;
      var cb = Math.sqrt(gn*gn - gmm) / denom;
      var ia = ic - n;      // (n-1,m)
      var ib = ia - n + 1;  // (n-2,m)
      p[ic] = ca * sinpp * p[ia] - cb * p[ib];
      q[ic] = ca * (sinpp * q[ia] + cospp * p[ia]) - cb * q[ib];
      ic++;
    }

    // sqrt(2n) p(n,n) = sqrt(2n-1) sqrt(1-x^2) p(n-1,n-1)
    var ca = Math.sqrt(1.0 - 0.5/n);
    var ia = ic - n - 1;  // (n-1,n-1)
    p[ic] = ca * cospp * p[ia];
    q[ic] = ca * (cospp * q[ia] - sinpp * p[ia]);
    ic++;
  }

  // Compute magnetic field vector components in geocentric frame.
  var x = 0;
  var y = 0;
  var z = 0;

  var ratio = 6371200.0 / r;
  var rr = ratio*ratio;

  ic = 1;
  for(var n = 1; n <= wmm_nmax; n++) {
    rr *= ratio;

    for(var m = 0; m <= n; m++) {
      var i = (ic-1) * 4;
    
      // Eq. 10-12.
      var g = (wmm_coef[i+0] + t * wmm_coef[i+2]);
      var h = (wmm_coef[i+1] + t * wmm_coef[i+3]);

      var cxz = (g * cl[m] + h * sl[m]) * rr;
      var cyy = (g * sl[m] - h * cl[m]) * rr;

      x -= cxz * q[ic];
      if(cospp != 0.0) {
        y += cyy * m * p[ic] / cospp;
      }
      else if(m == 1) {
        // Sect. 1.4: limits for p(n,m) / cospp as pp -> pi/2
        // m = 0  multiplies a zero coefficient
        // m = 1  finite limits
        // m > 1  -> 0
        // p(n,1) = cospp dp(n,0)/dx and
        // q(n,1) = -sinpp dp(n,0)/dx + cospp d^2p(n,0)/dx^2
        // so we can also write:
        // p(n,1) / cospp = -q(n,1) / sinpp when cospp -> 0
        y -= cyy * q[ic] / sinpp;
      }
      z -= (n+1.0) * cxz * p[ic];

      ic++;
    }
  }

  // Rotate into WGS84 frame.
  // x(WGS84) = x' cos(phi' - phi) - z' sin(phi' - phi)
  // y(WGS84) = y'
  // z(WGS84) = x' sin(phi' - phi) + z' cos(phi' - phi)
  var sindelt = sinpp * cosphi - cospp * sinphi;
  var cosdelt = cospp * cosphi + sinpp * sinphi;

  var result = [ x * cosdelt - z * sindelt,
                 y,
                 x * sindelt + z * cosdelt ];

  return result;
}
  

