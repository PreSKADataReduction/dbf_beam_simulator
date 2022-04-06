use std::f64::consts::PI;

use num::complex::Complex;

use crate::{
    constants::LIGHT_SPEED,
    utils::{calc_averaged_ant_output, integrate_az},
};

use scorus::{
    coordinates::{rotation3d::RotMatrix, SphCoord, Vec3d},
    healpix::{
        pix::pix2vec_ring,
        rotation::rotate_ring,
        utils::{npix2nside, nside2npix},
    },
};

#[allow(clippy::too_many_arguments)]
pub fn calc_array_beam(
    nside: usize,
    x_list: &[f64],
    y_list: &[f64],
    z_list: &[f64],
    w_list: &[f64],
    phi_list: &[f64],
    freq_Hz: f64,
    ground_cut: bool,
) -> Vec<f64> {
    let npix = nside2npix(nside);
    let lambda = LIGHT_SPEED / (freq_Hz);
    (0..npix)
        .map(|i| {
            if i < npix / 2 || !ground_cut {
                let pointing = pix2vec_ring::<f64>(nside, i);
                x_list
                    .iter()
                    .zip(
                        y_list
                            .iter()
                            .zip(z_list.iter().zip(w_list.iter().zip(phi_list.iter()))),
                    )
                    .map(|(&x, (&y, (&z, (&w, &phi))))| {
                        let dl = pointing[0] * x + pointing[1] * y + pointing[2] * z;
                        let phase = dl / lambda * 2.0 * PI;

                        Complex::from_polar(w, phase - phi)
                    })
                    .sum::<Complex<f64>>()
                    .norm_sqr()
            } else {
                0.0
            }
        })
        .collect()
}

pub fn calc_phase_from_pointing(
    x_list: &[f64],
    y_list: &[f64],
    z_list: &[f64],
    az_from_north: f64,
    zenith: f64,
    freq_Hz: f64,
) -> Vec<f64> {
    let az_from_x = PI / 2.0 - az_from_north;
    let lambda = LIGHT_SPEED / freq_Hz;
    let dir = Vec3d::from_sph_coord(SphCoord::new(zenith, az_from_x));
    x_list
        .iter()
        .zip(y_list.iter().zip(z_list.iter()))
        .map(|(&x, (&y, &z))| {
            let dl = dir[0] * x + dir[1] * y + dir[2] * z;
            dl / lambda * 2.0 * PI
        })
        .collect()
}

#[allow(clippy::too_many_arguments)]
pub fn calc_averaged_array_beam(
    lat_deg: f64,
    ant_beam: &[f64],
    x_list: &[f64],
    y_list: &[f64],
    z_list: &[f64],
    w_list: &[f64],
    phi_list: &[f64],
    freq_hz: f64,
) -> (Vec<f64>, Vec<usize>, Vec<f64>) {
    let nside = npix2nside(ant_beam.len());
    let array_beam = calc_array_beam(
        nside, x_list, y_list, z_list, w_list, phi_list, freq_hz, true,
    );
    let total_beam: Vec<f64> = array_beam
        .iter()
        .zip(ant_beam.iter())
        .map(|(&a, &b)| a * b)
        .collect();
    let rot =
        RotMatrix::about_axis_by_angle(&Vec3d::new(0.0, 1.0, 0.0), (90.0 - lat_deg).to_radians())
            * RotMatrix::about_axis_by_angle(&Vec3d::new(0.0, 0.0, 1.0), -90_f64.to_radians());
    let rotated_beam = rotate_ring(&total_beam, &rot);
    let (mean_beam, weight, theta) = integrate_az(&rotated_beam);
    (mean_beam, weight, theta)
}

#[allow(clippy::too_many_arguments)]
pub fn calc_averaged_ant_output2(
    sky: &[f64],
    lat_deg: f64,
    ant_beam: &[f64],
    x_list: &[f64],
    y_list: &[f64],
    z_list: &[f64],
    w_list: &[f64],
    phi_list: &[f64],
    freq_hz: f64,
) -> f64 {
    let (mean_beam, _weight, _theta) = calc_averaged_array_beam(
        lat_deg, ant_beam, x_list, y_list, z_list, w_list, phi_list, freq_hz,
    );
    calc_averaged_ant_output(&mean_beam, sky)
}
