// Base 10 to 60 conversion.

const unit_conv = {
  "rad": Math.PI / 180.0,
  "hr": 1.0 / 15.0,
  "deg": 1.0
};

function zeropad (value, width) {
  var tmpv = Math.abs(value).toString();
  var rpt = width - tmpv.length;

  var sign = "";
  if(value < 0) {
    sign = "-";
    rpt--;
  }
    
  if(rpt < 0)
    rpt = 0;
  
  var s = sign + "0".repeat(rpt) + tmpv;
  
  return s;
}

function spacepad (value, width) {
  var tmpv = value.toString();
  var rpt = width - tmpv.length;
  if(rpt < 0)
    rpt = 0;
  
  var s = " ".repeat(rpt) + tmpv;

  return s;
}

function base10_to_60(inangle, inunit, sep, sign, dp, outunit) {
  if(inangle < 0)
    sign = "-";

  inangle = Math.abs(inangle);

  var outangle = inangle * unit_conv[outunit] / unit_conv[inunit];

  var fs = 1;
  for(var i = 0; i < dp; i++)
    fs *= 10;

  var fm = fs * 60;
  var fh = fs * 3600;
  
  var ff = Math.round(outangle * fh);
  var hh = Math.floor(ff / fh);
  ff -= hh * fh;
  var mm = Math.floor(ff / fm);
  ff -= mm * fm;
  var ss = Math.floor(ff / fs);
  ff -= ss * fs;

  var result = sign + zeropad(hh, 2) + sep + zeropad(mm, 2) + sep + zeropad(ss, 2);
  if(dp > 0)
    result += "." + zeropad(ff, dp);

  return result;
}

// RA, Dec to standard coordinates.

function standc (a, d, tpa, tpd) {
  var sd = Math.sin(d);
  var cd = Math.cos(d);
  var stpd = Math.sin(tpd);
  var ctpd = Math.cos(tpd);
  
  var c = Math.cos(a - tpa);
  
  var denom = stpd * sd + ctpd * cd * c;
  
  var xi = Math.sin(a - tpa) * cd / denom;
  var xn = (ctpd * sd - stpd * cd * c) / denom;
  
  return [xi, xn];
}

// Standard coordinates to RA, Dec.

function xixn (xi, xn, tpa, tpd) {
  var stpd = Math.sin(tpd);
  var ctpd = Math.cos(tpd);
  
  var denom = ctpd - xn * stpd;
  
  var aa = Math.atan2(xi, denom);
  var a = (aa + tpa) % (2.0*Math.PI);
  if(a < 0)
    a += 2.0*Math.PI;
  
  var d = Math.atan2(stpd + xn * ctpd, Math.hypot(xi, denom));
  
  return [a, d];
}

// RA/DEC to unit vector.

function ad_to_v (a, d) {
  var sa = Math.sin(a);
  var ca = Math.cos(a);
  var sd = Math.sin(d);
  var cd = Math.cos(d);
  
  return [ ca*cd, sa*cd, sd ];
}

// Vector to RA/DEC.

function v_to_ad (n) {
  var a = Math.atan2(n[1], n[0]);
  var d = Math.atan2(n[2], Math.sqrt(n[0]*n[0]+n[1]*n[1]));

  if(a < 0)
    a += 2.0*Math.PI;
  
  return [ a, d ];
}

// N.B. these matrices are row-major 3x3 which is what I'm used to
// for positional astronomy.  They need to be converted for use
// with OpenGL.

function m_identity () {
  return [ 1.0, 0.0, 0.0,
           0.0, 1.0, 0.0,
           0.0, 0.0, 1.0 ];
}

function euler_rotate (m, axis, angle) {
  var sa = Math.sin(angle);
  var ca = Math.cos(angle);

  euler_rotate_sc(m, axis, sa, ca);
}

function euler_rotate_sc (m, axis, sa, ca) {
  var i, j;
  
  if(axis == 1) {
    i = 1;
    j = 2;
  }
  else if(axis == 2) {
    i = 2;
    j = 0;
  }
  else if(axis == 3) {
    i = 0;
    j = 1;
  }

  for(var k = 0; k < 3; k++) {
    var atmp =  ca * m[i*3+k] + sa * m[j*3+k];
    var btmp = -sa * m[i*3+k] + ca * m[j*3+k];
    m[i*3+k] = atmp;
    m[j*3+k] = btmp;
  }
}

function m_x_v (m, vi) {
  var vo = new Array(3);
  
  for(var j = 0; j < 3; j++) {
    vo[j] = 0;

    for(var i = 0; i < 3; i++)
      vo[j] += m[j*3+i] * vi[i];
  }

  return vo;
}

function mt_x_v (m, vi) {
  var vo = new Array(3);

  for(var j = 0; j < 3; j++) {
    vo[j] = 0;

    for(var i = 0; i < 3; i++)
      vo[j] += m[i*3+j] * vi[i];
  }

  return vo;
}

function m_x_m (a, b) {
  var c = new Array(9);

  for(var j = 0; j < 3; j++)
    for(var i = 0; i < 3; i++) {
      var tmp = 0;

      for(var k = 0; k < 3; k++)
        tmp += a[j*3+k] * b[k*3+i];

      c[j*3+i] = tmp;
    }

  return c;
}

function m_transpose (a) {
  var b = new Array(9);

  for(var j = 0; j < 3; j++)
    for(var i = 0; i < 3; i++)
      b[j*3+i] = a[i*3+j];

  return b;
}

// Convert 3x3 row-major to 4x4 column-major OpenGL matrices.

function m_to_gl (a) {
  var b = new Array(16);

  for(var j = 0; j < 3; j++) {
    for(var i = 0; i < 3; i++)
      b[j*4+i] = a[i*3+j];

    b[j*4+3] = 0;
  }

  for(var i = 0; i < 3; i++)
    b[3*4+i] = 0;

  b[3*4+3] = 1;

  return b;
}

// Precession and frame bias matrix, IAU 2006.  Time argument is TT Julian
// centuries since 2000.0.

function pfb_matrix (t) {
  // Precession and frame bias, IAU 2006, Fukushima-Williams angles.
  // Coefficients from Hilton et al. 2006.  Units are arcsec.
  var gam = ((((2.60e-8 * t - 2.788e-6) * t - 0.00031238) * t + 0.4932044) * t + 10.556378) * t - 0.052928;
  var phi = ((((-1.76e-8 * t - 4.40e-7) * t + 0.00053289) * t + 0.0511268) * t - 46.811016) * t + 84381.412819;
  var psi = ((((-1.48e-8 * t - 2.6452e-5) * t - 0.00018522) * t + 1.5584175) * t + 5038.481484) * t - 0.041775;
  var epsa = ((((-4.34e-8 * t - 5.76e-7) * t + 0.00200340) * t - 0.0001831) * t - 46.836769) * t + 84381.406000;

  var m = m_identity();
  euler_rotate(m, 3,  gam * Math.PI / (3600.0*180.0));
  euler_rotate(m, 1,  phi * Math.PI / (3600.0*180.0));
  euler_rotate(m, 3, -psi * Math.PI / (3600.0*180.0));
  euler_rotate(m, 1, -epsa * Math.PI / (3600.0*180.0));

  return m;
}

// Offset from Earth rotation angle to Greenwich Mean Sidereal Time,
// aka the slowly varying part.  IAU 2006, Capitaine et al. (2003).

function era_to_gmst (t) {
  var delt = ((((-0.0000000368 * t - 0.0000299560) * t - 0.00000044) * t + 1.39158170) * t + 4612.156534) * t + 0.014506;
  return(delt * Math.PI / (3600.0*180.0));
}

// Earth rotation angle, aka the rapidly varying part.
// Constants from Capitaine et al. 2000.

function era (jdtk) {
  var tmp = (0.7790572732640 + (jdtk % 1.0) + ((0.00273781191135448*jdtk) % 1.0)) % 1.0;
  return(tmp * 2.0 * Math.PI);
}

function lst_matrix (lst) {
  var lstm = m_identity();
  euler_rotate(lstm, 3, lst);
  return lstm;
}

function latitude_matrix (latitude) {
  var sinphi = Math.sin(latitude);
  var cosphi = Math.cos(latitude);

  var phm = m_identity();
  euler_rotate_sc(phm, 2, cosphi, sinphi);

  return phm;
}

// Get current MJD.

function getmjd () {
  var millis = Date.now();
  var dunix = millis * 1.0e-3 / 86400.0;
  return(dunix + 40587.0);
}

function mjd2jdtk (mjd) {
  return(mjd - 51544.5);
}

function jdtk2jctk (jdtk) {
  return(jdtk / (365.25 * 100.0));
}
