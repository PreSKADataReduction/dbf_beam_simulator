#![cfg(not(target_family = "wasm"))]
extern crate dbf_beam_simulator;

use clap::{Arg, ArgGroup, Command};

use ndarray::{s, Ix2};

use fitsimg::read_img;

use scorus::{
    coordinates::SphCoord,
    healpix::{
        interp::get_interpol_ring,
        utils::{npix2nside, nside2npix},
    },
};

use healpix_fits::{read_map, write_map};

use dbf_beam_simulator::regular_array::quarter_wgt2pattern;

fn main() {
    let matches = Command::new("wgt2beam")
        .arg(
            Arg::new("ant_beam")
                .short('a')
                .long("ant")
                .takes_value(true)
                .value_name("single antenna beam")
                .required(false)
                .help("ant beam"),
        )
        .arg(
            Arg::new("nside")
                .short('n')
                .long("nside")
                .takes_value(true)
                .value_name("nside")
                .required(false)
                .help("nside"),
        )
        .arg(
            Arg::new("wgt")
                .short('w')
                .long("wgt")
                .takes_value(true)
                .value_name("wgt file")
                .required(true)
                .help("wgt"),
        )
        .arg(
            Arg::new("spacing")
                .short('d')
                .long("spacing")
                .takes_value(true)
                .value_name("spacing in metre")
                .required(true)
                .help("spacing in metre"),
        )
        .arg(
            Arg::new("freq_MHz")
                .short('f')
                .long("freq")
                .takes_value(true)
                .value_name("freq in MHz")
                .required(true)
                .help("freq in MHz"),
        )
        .arg(
            Arg::new("out_beam")
                .short('o')
                .long("out")
                .takes_value(true)
                .value_name("output beam")
                .required(true)
                .help("output beam file"),
        )
        .group(
            ArgGroup::new("inputs")
                .args(&["nside", "ant_beam"])
                .required(true),
        )
        .get_matches();

    let (ant_beam, nside) = if let Some(fname) = matches.value_of("ant_beam") {
        let mut a = read_map::<f64>(fname, &["TEMPERATURE"], 1);
        let ant_beam = a.pop().unwrap();
        let nside = npix2nside(ant_beam.len());
        (ant_beam, nside)
    } else {
        let nside = matches.value_of("nside").unwrap().parse::<usize>().unwrap();
        let npix = nside2npix(nside);
        let ant_beam = (0..npix).map(|_| 1.0).collect();
        (ant_beam, nside)
    };

    let wgt = read_img::<f64>(matches.value_of("wgt").unwrap().to_string(), 0)
        .unwrap()
        .into_dimensionality::<Ix2>()
        .unwrap();
    let h = wgt.shape()[0];
    let w = wgt.shape()[1];

    let d = matches.value_of("spacing").unwrap().parse::<f64>().unwrap();

    let freq_mhz = matches
        .value_of("freq_MHz")
        .unwrap()
        .parse::<f64>()
        .unwrap();

    let array_beam =
        quarter_wgt2pattern(wgt.slice(s![h / 2..h, w / 2..w]).view(), d, freq_mhz, nside);

    let mut total_beam: Vec<_> = array_beam
        .iter()
        .zip(ant_beam.iter())
        .map(|(&a, &b)| a * b)
        .collect();

    let beam_norm = {
        let (p, w) = get_interpol_ring(nside, SphCoord::<f64>::new(0.0, 0.0));
        p.iter()
            .zip(w.iter())
            .map(|(&ipix, &w1)| total_beam[ipix] * w1)
            .sum::<f64>()
    };

    total_beam.iter_mut().for_each(|x| {
        *x /= beam_norm;
    });

    write_map(
        matches.value_of("out_beam").unwrap(),
        &[&total_beam],
        false,
        true,
    );
}
