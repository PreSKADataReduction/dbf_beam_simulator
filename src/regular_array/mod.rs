pub mod utils;
pub use utils::*;

use std::f64::consts::PI;

use num::complex::Complex;

use ndarray::{s, Array2, ArrayView2};

use crate::fft::{fft2, fftshift2};

use scorus::{
    coordinates::{SphCoord, Vec3d},
    healpix::{
        interp::get_interpol_ring,
        pix::pix2vec_ring,
        utils::{npix2nside, nside2npix},
    },
};

pub use crate::{
    arbitrary_array::{calc_array_beam, calc_averaged_array_beam, calc_phase_from_pointing},
    constants::LIGHT_SPEED,
    utils::{averaged_beam_to_healpix, integrate_az},
};

pub fn beam_opt_func_obj1(beam1: &[f64], beam2: &[f64], wgt: &[f64]) -> f64 {
    assert_eq!(beam1.len(), wgt.len());
    assert_eq!(beam2.len(), wgt.len());
    let npix = beam1.len();
    let domega = 4.0 * PI / npix as f64;
    let s1 = beam1
        .iter()
        .zip(wgt.iter())
        .map(|(&a, &b)| a * b)
        .sum::<f64>();
    let s2 = beam2
        .iter()
        .zip(wgt.iter())
        .map(|(&a, &b)| a * b)
        .sum::<f64>();
    beam1
        .iter()
        .zip(beam2.iter().zip(wgt.iter()))
        .map(|(&b1, (&b2, &w))| (b1 / s1 - b2 / s2) * w)
        .map(|x| x.powi(2))
        .sum::<f64>()
        / domega
}

#[allow(clippy::too_many_arguments)]
pub fn beam_opt_func_obj2(
    beam0: &[f64],
    lat_deg: f64,
    ant_beam: &[f64],
    x_list: &[f64],
    y_list: &[f64],
    z_list: &[f64],
    w_list: &[f64],
    phi_list: &[f64],
    freq_hz: f64,
) -> f64 {
    let (mean_beam, weight, _theta) = calc_averaged_array_beam(
        lat_deg, ant_beam, x_list, y_list, z_list, w_list, phi_list, freq_hz,
    );
    let weight: Vec<_> = weight.into_iter().map(|x| x as f64).collect();
    beam_opt_func_obj1(beam0, &mean_beam, &weight)
}

pub fn zenith_ns_sym_array(y: &[f64], w: &[f64]) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>) {
    let mut ant_x = vec![0.0];
    let mut ant_y = vec![0.0];
    let mut ant_z = vec![0.0];
    let mut ant_w = vec![1.0];

    for (&y1, &w1) in y.iter().zip(w.iter()) {
        ant_x.push(0.0);
        ant_y.push(y1);
        ant_z.push(0.0);
        ant_w.push(w1);

        ant_x.push(0.0);
        ant_y.push(-y1);
        ant_z.push(0.0);
        ant_w.push(w1);
    }
    (ant_x, ant_y, ant_z, ant_w)
}

pub fn sym_weight(w: &[f64]) -> Vec<f64> {
    let mut ant_w = vec![1.0];

    for &w1 in w {
        ant_w.push(w1);
        ant_w.push(w1);
    }
    ant_w
}

pub fn design_square_array(
    spacing: f64,
    freq_mega_hz: f64,
    sigma_deg: f64,
    n: isize,
) -> Array2<f64> {
    let center = n / 2;
    let lbd = LIGHT_SPEED / (freq_mega_hz * 1e6);
    let sigma = spacing / lbd * sigma_deg.to_radians().sin();
    let mut beam_pattern = Array2::<Complex<f64>>::zeros((n as usize, n as usize));
    for i in 0..n {
        let x = i - center;
        let fx = x as f64 / n as f64;
        let bx = f64::exp(-fx.powi(2) / (4.0 * sigma.powi(2))); //4 for sqrt(B)
        for j in 0..n {
            let y = j - center;
            let fy = y as f64 / n as f64;
            let by = f64::exp(-fy.powi(2) / (4.0 * sigma.powi(2))); //4 for sqrt(B)
            beam_pattern[(i as usize, j as usize)] = Complex::<f64>::from(bx * by);
        }
    }
    let mut beam_pattern = fftshift2(beam_pattern.view());
    let mut wgt = Array2::<Complex<f64>>::zeros((n as usize, n as usize));
    fft2(beam_pattern.view_mut(), wgt.view_mut());
    let norm = wgt[(0, 0)].re;
    let mut wgt = wgt.map(|x| x.re);
    wgt.iter_mut().for_each(|x| *x /= norm);
    println!("{}", wgt[(1, 1)]);
    println!("{}", wgt[(0, 0)]);
    fftshift2(wgt.view())
}

pub fn pattern2wgt(hp: &[f64], d: f64, freq_mhz: f64, n1: isize) -> Array2<f64> {
    let n = if n1 % 2 == 1 { n1 + 1 } else { n1 };
    let npix = hp.len();
    let nside = npix2nside(npix);
    let center = if n1 % 2 == 1 {
        n as f64 / 2.0
    } else {
        n as f64 / 2.0 - 0.5
    };
    let lbd = LIGHT_SPEED / (freq_mhz * 1e6);
    let u = d / lbd;

    let mut projected = Array2::<Complex<f64>>::zeros((n as usize, n as usize));

    for i in 0..n {
        let fx = (i as f64 - center) / n as f64;
        let nx = fx / u;
        for j in 0..n {
            let fy = (j as f64 - center) / n as f64;
            let ny = fy / u;
            let r2 = nx.powi(2) + ny.powi(2);
            if r2 > 1.0 {
                continue;
            }
            let nz = (1.0 - r2).sqrt();

            let (p, w) = get_interpol_ring(nside, SphCoord::from_xyz(nx, ny, nz));
            projected[(i as usize, j as usize)] = p
                .iter()
                .zip(w.iter())
                .map(|(&i1, &w1)| w1 * hp[i1])
                .sum::<f64>()
                .sqrt()
                .into();
        }
    }
    let mut projected = fftshift2(projected.view());
    let mut wgt = Array2::<Complex<f64>>::zeros((n as usize, n as usize));
    fft2(projected.view_mut(), wgt.view_mut());
    let norm = wgt[(0, 0)].re;
    let mut wgt = wgt.map(|x| x.re);
    wgt.iter_mut().for_each(|x| *x /= norm);
    let result = fftshift2(wgt.view());
    if n1 % 2 == 1 {
        result.slice(s![1.., 1..]).to_owned()
    } else {
        result
    }
}

pub fn wgt2pattern(wgt: ArrayView2<f64>, d: f64, freq_mhz: f64, nside: usize) -> Vec<f64> {
    let npix = nside2npix(nside);
    let lbd = LIGHT_SPEED / (freq_mhz * 1e6);
    let u = d / lbd;
    (0..npix)
        .map(|ipix| {
            if ipix < npix / 2 {
                let Vec3d { x: nx, y: ny, z: _ } = pix2vec_ring::<f64>(nside, ipix);
                let mut p = Complex::<f64>::new(0.0, 0.0);
                for i in 0..wgt.shape()[0] {
                    let m = i as f64 - (wgt.shape()[0] - 1) as f64 / 2.0;
                    for j in 0..wgt.shape()[1] {
                        let w = wgt[(i, j)];
                        let n = j as f64 - (wgt.shape()[1] - 1) as f64 / 2.0;
                        p += Complex::<f64>::from_polar(w, 2.0 * PI * (m * u * nx + n * u * ny));
                    }
                }
                p.norm_sqr()
            } else {
                0.0
            }
        })
        .collect()
}

pub fn quarter_wgt2pattern(
    quarter_wgt: ArrayView2<f64>,
    d: f64,
    freq_mhz: f64,
    nside: usize,
) -> Vec<f64> {
    let npix = nside2npix(nside);
    let lbd = LIGHT_SPEED / (freq_mhz * 1e6);
    let u = d / lbd;
    (0..npix)
        .map(|ipix| {
            if ipix < npix / 2 {
                let Vec3d { x: nx, y: ny, z: _ } = pix2vec_ring::<f64>(nside, ipix);
                let mut p = 0.0;
                for i in 0..quarter_wgt.shape()[0] {
                    for j in 0..quarter_wgt.shape()[1] {
                        let w = quarter_wgt[(i, j)];
                        p += w *  
                            (2.0 * PI * (i as f64) * u * nx).cos()
                                * (2.0 * PI * (j as f64) * u * ny).cos()
                                * if i == 0 { 1.0 } else { 2.0 }
                                * if j == 0 { 1.0 } else { 2.0 }
                        ;
                    }
                }
                p.powi(2)
            } else {
                0.0
            }
        })
        .collect()
}
