function make_ortho (xl, xh, yl, yh, zn, zf) {
  var dx = xh-xl;
  var dy = yh-yl;
  var dz = zf-zn;
    
  return [    2.0 / dx,         0.0,         0.0,   0.0,
                   0.0,    2.0 / dy,         0.0,   0.0,
                   0.0,         0.0,   -2.0 / dz,   0.0,
	   -(xh+xl)/dx, -(yh+yl)/dy, -(zf+zn)/dz,   1.0  ];
}

function make_frustrum (xl, xh, yl, yh, zn, zf) {
  var dx = xh-xl;
  var dy = yh-yl;
  var dz = zf-zn;

  return [   2.0*zn/dx,        0.0,           0.0,  0.0,
	           0.0,  2.0*zn/dy,           0.0,  0.0,
            (xh+xl)/dx, (yh+yl)/dy,   -(zf+zn)/dz, -1.0,
                   0.0,        0.0, -2.0*zf*zn/dz,  0.0  ];
}

// tangent is tan of half the vertical angle of view.

function make_perspective (tangent, aspect, zn, zf) {
  var height = zn * tangent;
  var width = height * aspect;

  return make_frustrum(-width, width, -height, height, zn, zf);
}

function make_ortho_fov (tangent, aspect, zn, zf) {
  var height = zn * tangent;
  var width = height * aspect;

  return make_ortho(-width, width, -height, height, zn, zf);
}

function renorm (v) {
  var mod = Math.sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
  if(mod > 0) {
    v[0] /= mod;
    v[1] /= mod;
    v[2] /= mod;
  }
}

function cross_product (a, b) {
  var c = new Array(3);

  c[0] = a[1]*b[2] - a[2]*b[1];
  c[1] = a[2]*b[0] - a[0]*b[2];
  c[2] = a[0]*b[1] - a[1]*b[0];

  return c;
}

function look_at (cam, targ, up) {
  var plusz = new Array(3);

  // Unit vector from camera to target (-z).
  plusz[0] = cam[0] - targ[0];
  plusz[1] = cam[1] - targ[1];
  plusz[2] = cam[2] - targ[2];
  renorm(plusz);

  // Unit vector pointing +x in camera view.
  var plusx = cross_product(up, plusz);
  renorm(plusx);

  // Unit vector pointing +y in camera view.
  var plusy = cross_product(plusz, plusx);

  // Build "look at" matrix.
  return [ plusx[0], plusy[0], plusz[0], 0.0,
           plusx[1], plusy[1], plusz[1], 0.0,
           plusx[2], plusy[2], plusz[2], 0.0,
          -plusx[0]*cam[0]-plusx[1]*cam[1]-plusx[2]*cam[2], // -x . cam
          -plusy[0]*cam[0]-plusy[1]*cam[1]-plusy[2]*cam[2], // -y . cam
          -plusz[0]*cam[0]-plusz[1]*cam[1]-plusz[2]*cam[2], // -z . cam
           1.0 ];
}

function matrix_multiply (a, b) {
  var c = new Array(16);

  for(var i = 0; i < 4; i++) {
    for(var j = 0; j < 4; j++) {
      var sum = 0;

      for(k = 0; k < 4; k++)
        sum += a[k*4+j] * b[i*4+k];

      c[i*4+j] = sum;
    }
  }

  return c;
}

function matrix_transpose (a) {
  var b = new Array(16);

  for(var i = 0; i < 4; i++) {
    for(var j = 0; j < 4; j++) {
      b[i*4+j] = a[j*4+i];
    }
  }

  return b;
}

function forward_transform (a, b) {
  var c = new Array(4);

  for(var j = 0; j < 4; j++) {
    var sum = 0;

    for(k = 0; k < 4; k++)
      sum += a[k*4+j] * b[k];

    c[j] = sum;
  }

  return c;
}

function reverse_transform (a, b) {
  var c = new Array(4);

  for(var j = 0; j < 4; j++) {
    var sum = 0;

    for(k = 0; k < 4; k++)
      sum += a[j*4+k] * b[k];

    c[j] = sum;
  }

  return c;
}

function load_shader (gl, type, source) {
  const shader = gl.createShader(type);
  gl.shaderSource(shader, source);
  gl.compileShader(shader);

  if(!gl.getShaderParameter(shader, gl.COMPILE_STATUS)) {
    console.log("Shader compilation error: " + gl.getShaderInfoLog(shader) + source);
    gl.deleteShader(shader);
    return null;
  }

  return shader;
}

function load_program (gl, src_vertex_shader, src_fragment_shader) {
  const vertex_shader = load_shader(gl, gl.VERTEX_SHADER, src_vertex_shader);
  if(vertex_shader === null)
    return null;

  const fragment_shader = load_shader(gl, gl.FRAGMENT_SHADER, src_fragment_shader);
  if(fragment_shader === null) {
    glDeleteShader(vertex_shader);
    return null;
  }

  const shader_program = gl.createProgram();
  gl.attachShader(shader_program, vertex_shader);
  gl.attachShader(shader_program, fragment_shader);
  gl.linkProgram(shader_program);

  if(!gl.getProgramParameter(shader_program, gl.LINK_STATUS)) {
    console.log("Unable to initialize the shader program: " + gl.getProgramInfoLog(shader_program));
    glDeleteProgram(shader_program);
    glDeleteShader(vertex_shader);
    glDeleteShader(fragment_shader);
    return null;
  }

  var sh = { program: shader_program };

  const nattr = gl.getProgramParameter(shader_program, gl.ACTIVE_ATTRIBUTES);
  for(var iattr = 0; iattr < nattr; iattr++) {
    const info = gl.getActiveAttrib(shader_program, iattr);
    const loc = gl.getAttribLocation(shader_program, info.name);
    sh[info.name] = loc;
  }

  const nunif = gl.getProgramParameter(shader_program, gl.ACTIVE_UNIFORMS);
  for(var iunif = 0; iunif < nunif; iunif++) {
    const info = gl.getActiveUniform(shader_program, iunif);

    if(info.size > 1) {
      var pos = info.name.indexOf("[");
      var basename = info.name.substr(0, pos);

      var ent = new Array(info.size);

      for(var i = 0; i < info.size; i++) {
        ent[i] = gl.getUniformLocation(shader_program, basename+"["+i+"]");
      }

      sh[basename] = ent;
    }
    else {
      const loc = gl.getUniformLocation(shader_program, info.name);
      sh[info.name] = loc;
    }
  }

  return sh;
}

function makeellipsoid (xcen, ycen, zcen, rx, ry, rz, ni, nj, vertices) {
  var di = Math.PI / ni;
  var dj = 2*Math.PI / (nj-1);

  for(var i = 0; i < ni; i++) {
    var ai = i*di;    
    for(var j = 0; j < nj; j++) {
      var aj = j*dj;

      var x = rx * Math.cos(aj) * Math.sin(ai);
      var y = ry * Math.sin(aj) * Math.sin(ai);
      var z = rz * Math.cos(ai);

      var norm = 1.0/Math.sqrt(x*x+y*y+z*z);

      var ioff = (i*nj+j)*12;
      
      vertices[ioff+0] = xcen+x;
      vertices[ioff+1] = ycen+y;
      vertices[ioff+2] = zcen+z;
      vertices[ioff+3] = x*norm;
      vertices[ioff+4] = y*norm;
      vertices[ioff+5] = z*norm;

      var x = rx * Math.cos(aj) * Math.sin(ai+di);
      var y = ry * Math.sin(aj) * Math.sin(ai+di);
      var z = rz * Math.cos(ai+di);

      vertices[ioff+6] = xcen+x;
      vertices[ioff+7] = ycen+y;
      vertices[ioff+8] = zcen+z;
      vertices[ioff+9] = x*norm;
      vertices[ioff+10] = y*norm;
      vertices[ioff+11] = z*norm;
    }
  }
}
