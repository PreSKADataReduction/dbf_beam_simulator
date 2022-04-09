#![cfg(not(target_family = "wasm"))]
extern crate dbf_beam_simulator;

use clap::{Arg, Command};

use ndarray::s;

use fitsimg::write_img;

use scorus::{
    coordinates::SphCoord,
    healpix::{interp::get_interpol_ring, utils::npix2nside},
};

use healpix_fits::read_map;

use dbf_beam_simulator::regular_array::pattern2wgt;

fn main() {
    let matches = Command::new("gaussian beam")
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
            Arg::new("target_beam")
                .short('t')
                .long("tb")
                .takes_value(true)
                .value_name("target beam")
                .required(true)
                .help("target_beam"),
        )
        .arg(
            Arg::new("array_size")
                .short('y')
                .long("as")
                .takes_value(true)
                .value_name("array_size")
                .required(true)
                .help("array size"),
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
            Arg::new("out_wgt")
                .short('o')
                .long("out")
                .takes_value(true)
                .value_name("output wgt file name")
                .required(true)
                .help("output file name"),
        )
        .get_matches();

    let target_beam = read_map::<f64>(
        matches.value_of("target_beam").unwrap(),
        &["TEMPERATURE"],
        1,
    )
    .pop()
    .unwrap();

    let npix = target_beam.len();
    let nside = npix2nside(npix);

    let ant_beam = if let Some(fname) = matches.value_of("ant_beam") {
        let ant_beam = read_map::<f64>(fname, &["TEMPERATURE"], 1).pop().unwrap();
        assert_eq!(npix, ant_beam.len());
        ant_beam
    } else {
        (0..npix).map(|_| 1.0).collect()
    };

    let array_size = matches
        .value_of("array_size")
        .unwrap()
        .parse::<isize>()
        .unwrap();
    let d = matches.value_of("spacing").unwrap().parse::<f64>().unwrap();

    let freq_mhz = matches
        .value_of("freq_MHz")
        .unwrap()
        .parse::<f64>()
        .unwrap();

    let mut array_beam: Vec<_> = ant_beam
        .iter()
        .zip(target_beam.iter())
        .map(|(&a, &t)| t / a)
        .collect();

    let beam_norm = {
        let (p, w) = get_interpol_ring(nside, SphCoord::<f64>::new(0.0, 0.0));
        p.iter()
            .zip(w.iter())
            .map(|(&ipix, &w1)| array_beam[ipix] * w1)
            .sum::<f64>()
    };

    array_beam.iter_mut().for_each(|x| {
        *x /= beam_norm;
        if *x > 2.0 {
            eprintln!("warning beam max >2: {}", *x);
        }
    });

    let mut wgt = pattern2wgt(&array_beam, d, freq_mhz, array_size);
    wgt.slice_mut(s![0, ..]).fill(0.0);
    wgt.slice_mut(s![.., 0]).fill(0.0);
    //let pattern2=quarter_wgt2pattern(
    //     wgt.slice(s![array_size/2..array_size, array_size/2..array_size]).view()
    //     , d, freq_mhz, nside);
    write_img(
        matches.value_of("out_wgt").unwrap().to_string(),
        &wgt.into_dyn(),
    )
    .unwrap();
}
