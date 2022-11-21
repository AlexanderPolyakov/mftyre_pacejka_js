var trace1 = {
  x: [1, 2, 3, 4],
  y: [10, 15, 13, 17],
  type: 'scatter'
};

function sin_mf(D, C, B, E) {
  return D * Math.sin(C * Math.atan(B - E * (B - Math.atan(B))))
}

function cos_mf(D, C, B, E) {
  return D * Math.cos(C * Math.atan(B - E * (B - Math.atan(B))))
}

// table A3.1
const p = {
  // longitudinal
  Cx1: 1.579,
  Dx1: 1.0422,
  Dx2: -0.08285,
  Dx3: 0,
  Ex1: 0.11113,
  Ex2: 0.3143,
  Ex3: 0,
  Ex4: 0.001719,
  Hx1: 2.1615e-4,
  Hx2: 0.0011598,
  Kx1: 21.687,
  Kx2: 13.728,
  Kx3: -0.4098,
  Vx1: 2.0283e-5,
  Vx2: 1.0568e-4,
  px1: -0.3485,
  px2: 0.37824,
  px3: -0.09603,
  px4: 0.06518,

  // lateral
  Cy1: 1.338,
  Dy1: 0.8785,
  Dy2: -0.06452,
  Dy3: 0,
  Ey1: -0.8057,
  Ey2: -0.6046,
  Ey3: 0.09854,
  Ey4: -6.697,
  Ey5: 0,
  Ky1: -15.324,
  Ky2: 1.715,
  Ky3: 0.3695,
  Ky4: 2.0005,
  Ky5: 0,
  Hy1: -0.001806,
  Hy2: 0.00352,
  Vy1: -0.00661,
  Vy2: 0.03592,
  Vy3: -0.162,
  Vy4: -0.4864,
  Ky6: -0.8987,
  Ky7: -0.4864,
  py1: -0.6255,
  py2: -0.06523,
  py3: -0.16666,
  py4: 0.2811,
  py5: 0
};

const lambda = {
  // longitudinal
  Cx: 1,
  Ex: 1,
  Kxk: 1,
  mux: 1,
  muV: 0, // set to 0 as 'not used'
  Hx: 1,
  Vx: 1,
  Fz0: 1,

  // lateral
  Cy: 1,
  Ey: 1,
  muy: 1,
  Hy: 1,
  Vy: 1,
  K_yalpha: 1,
  K_ygamma: 1,
};

// TODO: move into input
const ksi_0 = 1.0;
const ksi_1 = 1.0;
const ksi_2 = 1.0;
const ksi_3 = 1.0;
const ksi_4 = 1.0;

function calc_df_z(params, input) {
  const F_z0_prime = lambda.Fz0 * params.F_z0; // 4.E1
  const df_z = (input.F_z - F_z0_prime) / F_z0_prime; // 4.E2a
  return df_z;
}

function calc_dp_i(params, input) {
  const dp_i = (input.p_i - params.p_i0) / params.p_i0; // 4.E2b
  return dp_i;
}

function calc_mu_star(lambda_mu, params, input) {
  const lambda_mu_star = lambda_mu / (1 + lambda.muV * input.V_s/params.V_0); // 4.E7
  return lambda_mu_star;
}

function calc_mu_prime(lambda_mu_star) {
  const A_mu = 10;
  const lambda_mu_prime = A_mu * lambda_mu_star / (1 + (A_mu - 1) * lambda_mu_star); // 4.E8
  return lambda_mu_prime;
}

function longitudinal(params, input) {
  const eps_x = 0.1;

  const df_z = calc_df_z(params, input);
  const dp_i = calc_dp_i(params, input);

  const lambda_mux_star = calc_mu_star(lambda.mux, params, input);
  const lambda_mux_prime = calc_mu_prime(lambda_mux_star);

  const S_Hx = (p.Hx1 + p.Hx2 * df_z) * lambda.Hx; // 4.E17
  const kappa_x = input.kappa + S_Hx; // 4.E10
  const C_x = p.Cx1 * lambda.Cx; // 4.E11
  const mu_x = (p.Dx1 + p.Dx2 * df_z) * (1 + p.px3 * dp_i + p.px4 * dp_i ** 2) * (1 - p.Dx3 * input.gamma ** 2) * lambda_mux_star; // 4.E13
  const E_x = (p.Ex1 + p.Ex2 * df_z + p.Ex3 * df_z ** 2) * (1 - p.Ex4 * Math.sign(kappa_x)) * lambda.Ex; // 4.E14
  const K_xk = input.F_z * (p.Kx1 + p.Kx2 * df_z) * Math.exp(p.Kx3 * df_z) * (1 + p.px1 * dp_i + p.px2 * dp_i ** 2) * lambda.Kxk; // 4.E15
  const D_x = mu_x * input.F_z * ksi_1; // 4.E12
  const B_x = K_xk / (C_x * D_x + eps_x); // 4.E16
  const S_Vx = input.F_z * (p.Vx1 + p.Vx2 * df_z) * lambda.Vx * lambda_mux_prime * ksi_1; // 4.E18
  const F_x0 = sin_mf(D_x, C_x, B_x * kappa_x, E_x) + S_Vx; // 4.E9
 
  return F_x0;
}

function lateral(params, input) {
  const eps_y = 0.1;
  const eps_K = 0.1;

  const df_z = calc_df_z(params, input);
  const dp_i = calc_dp_i(params, input);

  const lambda_muy_star = calc_mu_star(lambda.muy, params, input);
  const lambda_muy_prime = calc_mu_prime(lambda_muy_star);

  const F_z0_prime = lambda.Fz0 * params.F_z0; // 4.E1
  const alpha_star = Math.tan(input.alpha) * Math.sign(input.V_cx); // 4.E3
  const gamma_star = Math.sin(input.gamma); // 4.E4

  const K_yalpha = p.Ky1 * F_z0_prime * (1 + p.py1 * dp_i) * (1 - p.Ky3 * Math.abs(gamma_star)) *
    Math.sin(p.Ky4 * Math.atan(input.F_z / F_z0_prime / ((p.Ky2 + p.Ky5 * gamma_star ** 2) * (1 + p.py2 * dp_i)))) * ksi_3 * lambda.K_yalpha; // 4.E25
  const C_y = p.Cy1 * lambda.Cy; // 4.E21
  const mu_y = (p.Dy1 + p.Dy2 * df_z) * (1 + p.py3 * dp_i + p.py4 * dp_i ** 2) * (1 - p.Dy3 * gamma_star ** 2) * lambda_muy_star; // 4.E23
  const D_y = mu_y * input.F_z * ksi_2; // 4.E22
  const B_y = K_yalpha / (C_y * D_y + eps_y); // 4.E26
  const S_Vygamma = input.F_z * (p.Vy3 + p.Vy4 * df_z) * gamma_star * lambda.K_ygamma * lambda_muy_prime * ksi_2; // 4.E28
  const S_Vy = input.F_z * (p.Vy1 + p.Vy2 * df_z) * lambda.Vy * lambda_muy_prime * ksi_2 + S_Vygamma; // 4.E29
  const K_ygamma0 = input.F_z * (p.Ky6 + p.Ky7 * df_z) * (1 + p.py5 * dp_i) * lambda.K_ygamma; // 4.E30
  const S_Hy = (p.Hy1 + p.Hy2 * df_z) * lambda.Hy + (K_ygamma0 * gamma_star - S_Vygamma) / (K_yalpha + eps_K) * ksi_0 + ksi_4 - 1; // 4.E27
  const alpha_y = alpha_star + S_Hy; // 4.E20
  const E_y = (p.Ey1 + p.Ey2 * df_z) * (1 + p.Ey5 * gamma_star ** 2 - (p.Ey3 + p.Ey4 * gamma_star) * Math.sign(alpha_y)) * lambda.Ey; // 4.E24
  const F_y0 = sin_mf(D_y, C_y, B_y * alpha_y, E_y) + S_Vy; // 4.E19

  return F_y0;
}

let longitudinal_f = {
  x: [],
  y: [],
  type: 'scatter',
  name: 'Fx0'
};

const defParams = {
  p_i0: 2.17,
  V_0: 16.67,
  F_z0: 4000
}

for (let kappa = -1; kappa < 1; kappa += 0.01) {
  longitudinal_f.x.push(kappa);
  longitudinal_f.y.push(longitudinal(defParams,
    {
      F_z: 4000,
      p_i: 2.17,
      kappa: kappa,
      V_s: 16.67,
      gamma: 0
    }));
}

let lateral_f = {
  x: [],
  y: [],
  type: 'scatter',
  name: 'Fy0'
};


for (let alpha = -90; alpha < 91; alpha += 0.5) {
  const alphaRad = alpha / 180 * Math.PI;
  lateral_f.x.push(alpha);
  lateral_f.y.push(lateral(defParams,
    {
      F_z: 4000,
      p_i: 2.17,
      kappa: 0,
      V_s: 16.67,
      V_cx: 16.67,
      gamma: 0,
      alpha: alphaRad
    }));
}

var data_long_f = [longitudinal_f];

Plotly.newPlot('fx0_plot', data_long_f,
  {
    title: "longitudinal",
    xaxis: {
      title: 'kappa'
    },
    yaxis: {
      title: 'Fx0'
    }
  });

Plotly.newPlot('fy0_plot', [lateral_f],
  {
    title: "lateral",
    xaxis: {
      title: 'alpha'
    },
    yaxis: {
      title: 'Fy0'
    }
  });

