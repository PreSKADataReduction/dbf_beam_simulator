use scorus::healpix::{
    pix::{pix2ring_ring, ring2z_ring},
    utils::{npix2nside, nring2nside, nside2npix, nside2nring},
};

pub fn integrate_az(hmap: &[f64]) -> (Vec<f64>, Vec<usize>, Vec<f64>) {
    let npix = hmap.len();
    let nside = npix2nside(npix);
    let nring = nside2nring(nside);
    let mut wgt = vec![0_usize; nring];
    let mut mean_values = vec![0.0; nring];
    let ring_idx: Vec<_> = (0..npix).map(|ipix| pix2ring_ring(nside, ipix)).collect();
    ring_idx.iter().for_each(|&i| wgt[i - 1] += 1);
    hmap.iter().zip(ring_idx.iter()).for_each(|(&x, &iring)| {
        mean_values[iring - 1] += x;
    });
    mean_values
        .iter_mut()
        .zip(wgt.iter())
        .for_each(|(x, &w)| *x /= w as f64);

    let theta = (1..=nring)
        .map(|iring| ring2z_ring::<f64>(nside, iring).acos())
        .collect();
    (mean_values, wgt, theta)
}

pub fn averaged_beam_to_healpix(beam: &[f64]) -> Vec<f64> {
    let nring = beam.len();
    let nside = nring2nside(nring);
    let npix = nside2npix(nside);
    (0..npix)
        .map(|ipix| beam[pix2ring_ring(nside, ipix) - 1])
        .collect()
}

pub fn calc_averaged_ant_output(beam: &[f64], sky: &[f64]) -> f64 {
    let beam_hp = averaged_beam_to_healpix(beam);
    assert_eq!(beam_hp.len(), sky.len());
    let norm = beam_hp.iter().sum::<f64>();
    beam_hp
        .iter()
        .zip(sky.iter())
        .map(|(&b, &s)| b * s)
        .sum::<f64>()
        / norm
}

pub use crate::arbitrary_array::calc_averaged_ant_output2;
