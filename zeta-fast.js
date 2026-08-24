function cgammaphi(x, y) { const z = toComplex(x, y); return cexp(cmul(z, complex(-3.53102424697, 1.57079632679).add(clog(complex(1.0, 0.0).add(z.scale(-2.0))))).add(clog(z).scale(-0.5)).add(complex(1.4189385332, 0.0))); }

function czeta_helper_4(z, bases, ns) {
  const z_ = toComplex(z);
  const cutoff = sqrt(z_.y * 0.159154943092) + 100.0 * step(z_.y, 120.0);
  const cutoffscale = 0.884 * cutoff + 0.442;
  let resX = 0.0, resY = 0.0, resX2 = 0.0, resY2 = 0.0;
  for (let row = 0; row < 4; row++) {
    const mags = [];
    const mags2 = [];
    const re = [];
    const im = [];
    for (let col = 0; col < 4; col++) {
      const mag = Math.exp(z_.x * bases[row][col]);
      const cutoffv = clamp(cutoffscale - 0.884 * ns[row][col], 0.0, 1.0);
      mags.push(mag * cutoffv);
      mags2.push(cutoffv / (ns[row][col] * mag));
      re.push(Math.cos(z_.y * bases[row][col]));
      im.push(Math.sin(z_.y * bases[row][col]));
    }
    let dotRe = 0.0, dotIm = 0.0, dotRe2 = 0.0, dotIm2 = 0.0;
    for (let col = 0; col < 4; col++) {
      dotRe += mags[col] * re[col];
      dotIm += mags[col] * im[col];
      dotRe2 += mags2[col] * re[col];
      dotIm2 += mags2[col] * im[col];
    }
    resX += dotRe;
    resY += dotIm;
    resX2 += dotRe2;
    resY2 += -dotIm2;
  }
  return { x: resX, y: resY, x2: resX2, y2: resY2 };
}

function czeta_helper_2(z, bases, coeffs) {
  const z_ = toComplex(z);  
  let resX = 0.0, resY = 0.0;
  for (let row = 0; row < 4; row++) {
    const mags = [];
    const re = [];
    const im = [];
    for (let col = 0; col < 4; col++) {
      mags.push(Math.exp(z_.x * bases[row][col]) * coeffs[row][col]);
      re.push(Math.cos(z_.y * bases[row][col]));
      im.push(Math.sin(z_.y * bases[row][col]));
    }
    for (let col = 0; col < 4; col++) {
      resX += mags[col] * re[col];
      resY += mags[col] * im[col];
    }
  }
  return complex(resX, resY);
}

function czeta_strip(z) {
  const z_ = toComplex(z);
  const bases1 = [
    [-0.69314718056, -1.09861228866811, -1.38629436111989, -1.6094379124341],
    [-1.79175946922805, -1.94591014905531, -2.07944154167984, -2.19722457733622],
    [-2.30258509299405, -2.39789527279837, -2.48490664978800, -2.56494935746154],
    [-2.63905732961526, -2.70805020110221, -2.77258872223978, -2.83321334405622]
  ];
  const ns1 = [
    [2.0, 3.0, 4.0, 5.0],
    [6.0, 7.0, 8.0, 9.0],
    [10.0, 11.0, 12.0, 13.0],
    [14.0, 15.0, 16.0, 17.0]
  ];
  const bases2 = [
    [-2.89037175789616, -2.94443897916644, -2.99573227355399, -3.04452243772342],
    [-3.09104245335832, -3.13549421592915, -3.17805383034795, -3.21887582486820],
    [-3.25809653802148, -3.29583686600433, -3.33220451017520, -3.36729582998647],
    [-3.40119738166216, -3.43398720448515, -3.46573590279973, -3.49650756146648]
  ];
  const ns2 = [
    [18.0, 19.0, 20.0, 21.0],
    [22.0, 23.0, 24.0, 25.0],
    [26.0, 27.0, 28.0, 29.0],
    [30.0, 31.0, 32.0, 33.0]
  ];
  const zeta_est1 = czeta_helper_4(z_, bases1, ns1);
  const zeta_est2 = czeta_helper_4(z_, bases2, ns2);
  const zetaX = 1.0 + zeta_est1.x + zeta_est2.x;
  const zetaY = 0.0 + zeta_est1.y + zeta_est2.y;
  const zetaX2 = 1.0 + zeta_est1.x2 + zeta_est2.x2;
  const zetaY2 = 0.0 + zeta_est1.y2 + zeta_est2.y2;
  const zetaA = complex(zetaX, zetaY);
  const zetaB = cdot(complex(zetaX2, zetaY2), cgammaphi(complex(1.0 - z_.x, z_.y)));
  if (z_.y < 120.0) {
    const t = 1.0 - Math.min(z_.x, 1.0);
    const alpha = t * t * (3.0 - 2.0 * t);
    return cadd(cscale(zetaA, 1.0 - alpha), cscale(zetaB, alpha));
  }
  return cadd(zetaA, zetaB);
}

function ceta_right(z) {
  const z_ = toComplex(z);
  const one = complex(1.0, 0.0);
  if (z_.x < 3.0 && z_.y > 54.0) {
    return cmul(csub(one, cscale(cexp(cscale(z_, -0.69314718056)), 2.0)), czeta_strip(z_));
  }
  let result = one;
  const bases1 = [
    [-0.69314718056, -1.09861228866811, -1.38629436111989, -1.6094379124341],
    [-1.79175946922805, -1.94591014905531, -2.07944154167984, -2.19722457733622],
    [-2.30258509299405, -2.39789527279837, -2.48490664978800, -2.56494935746154],
    [-2.63905732961526, -2.70805020110221, -2.77258872223978, -2.83321334405622]
  ];
  const coeffs1 = [
    [-1.00000000000000, 1.00000000000000, -1.00000000000000, 1.00000000000000],
    [-0.99999999999995, 0.99999999999847, -0.99999999996425, 0.99999999937104],
    [-0.99999999142280, 0.99999990708781, -0.99999918494666, 0.99999411949279],
    [-0.99996466193028, 0.99982127062071, -0.99923254216349, 0.99718148818347]
  ];
  const bases2 = [
    [-2.89037175789616, -2.94443897916644, -2.99573227355399, -3.04452243772342],
    [-3.09104245335832, -3.13549421592915, -3.17805383034795, -3.21887582486820],
    [-3.25809653802148, -3.29583686600433, -3.33220451017520, -3.36729582998647],
    [-3.40119738166216, -3.43398720448515, -3.46573590279973, -3.49650756146648]
  ];
  const coeffs2 = [
    [-0.99109047939434, 0.97562125072353, -0.94195422388664, 0.87910910712444],
    [-0.77852772396727, 0.64073335549827, -0.47964042231228, 0.31968999219853],
    [-0.18572334624203, 0.09196690020310, -0.03784892366350, 0.01254701255408],
    [-0.00320994916827, 0.00059346134942, -0.00007044051625, 0.00000402517236]
  ];
  result = cadd(result, czeta_helper_2(z_, bases1, coeffs1));
  result = cadd(result, czeta_helper_2(z_, bases2, coeffs2));
  return result;
}

function ceta_left(w) {
  const z_ = toComplex(-w.x, w.y);
  let component_a;
  const log_r = Math.log(z_.x * z_.x + z_.y * z_.y) * 0.5;
  if (z_.y > 200.0) {
    const theta = 1.57079632679;
    component_a = cmul_i(cexp(complex(z_.x + (log_r - 1.0) * z_.x - 0.5 * log_r, (log_r - 1.0) * z_.y - 0.5 * theta))).scale(1.2533141373155);
  } else if (z_.y > 20.0) {
    const theta = Math.atan2(z_.y, z_.x);
    component_a = cmul_i(cexp(cmul_i(z_).scale(theta - 1.57079632679).add(complex((log_r - 1.0) * z_.x - 0.5 * log_r, (log_r - 1.0) * z_.y - 0.5 * theta)))).scale(1.2533141373155);
  } else {
    component_a = cmul(cgamma(z_), csin(z_.scale(1.57079632679)));
  }
  const zadd1 = complex(z_.x + 1.0, z_.y);
  const two_neg_z = cexp(-0.69314718056 * z_);
  const multiplier = cmul(cexp(zadd1.scale(-1.14472988585)), cmul(complex(1.0 - 0.5 * two_neg_z.x, -0.5 * two_neg_z.y), cinv(complex(1.0 - two_neg_z.x, -two_neg_z.y))));
  const component_b = cmul(z_, ceta_right(zadd1));
  return 2.0 * conj(cmul(cmul(component_a, component_b), multiplier));
}

function ceta(x, y) {
  let z = toComplex(x, y);
  let conjugate_mask = complex(1.0, 1.0);
  if (z.y < 0.0) { z = complex(z.x, -z.y); conjugate_mask = complex(1.0, -1.0); }
  return cmulcomp(z.x < 0.0 ? ceta_left(z) : ceta_right(z), conjugate_mask);
}

function czeta(x, y) { const z = toComplex(x, y); return cmul(ceta(z), cinv(cadd(complex(1.0, 0.0), cexp(z.scale(-0.69314718056)).scale(-2.0)))); }
