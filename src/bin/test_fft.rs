#![cfg(not(target_family = "wasm"))]
#![allow(unused_imports)]
extern crate dbf_beam_simulator;

use ndarray::s;

use fitsimg::write_img;

use scorus::healpix::{pix::pix2ang_ring, utils::nside2npix};

use healpix_fits::write_map;

use dbf_beam_simulator::regular_array::{pattern2wgt, quarter_wgt2pattern, wgt2pattern};

fn main() {
    /*
    let nside = 128;
    let array_size = 8;
    let d = 1.5;
    let freq_mhz = 80.0;
    let npix = nside2npix(nside);
    let sigma: f64 = 15.0;
    let hpmap: Vec<_> = (0..npix)
        .map(|i| {
            (-pix2ang_ring::<f64>(nside, i).pol.to_degrees().powi(2) / (2.0 * sigma.powi(2))).exp()
        })
        .collect();
    write_map("input_pattern.fits", &[&hpmap], false, true);
    let mut wgt = pattern2wgt(&hpmap, d, freq_mhz, array_size);
    wgt.slice_mut(s![0, ..]).fill(0.0);
    wgt.slice_mut(s![.., 0]).fill(0.0);

    let pattern1 = wgt2pattern(wgt.view(), d, freq_mhz, nside);
    let pattern2 = quarter_wgt2pattern(
        wgt.slice(s![array_size / 2..array_size, array_size / 2..array_size])
            .view(),
        d,
        freq_mhz,
        nside,
    );
    write_img("wgt.fits".to_string(), &wgt.into_dyn()).unwrap();
    write_map("pattern1.fits", &[&pattern1], false, true);
    write_map("pattern2.fits", &[&pattern2], false, true);
    */
}
