extern crate dbf_beam_simulator;

use clap::{
    App
    ,Arg
};


use scorus::{
    healpix::{
        utils::{
            nside2npix
        }
        , pix::{
            pix2ang_ring
        }
    }
};

use healpix_fits::{
    write_map
};

fn main(){
    let matches=App::new("gaussian beam")
    .arg(
        Arg::new("nside")
        .short('n')
        .long("nside")
        .takes_value(true)
        .value_name("nside")
        .required(true)
        .help("nside")
    )
    .arg(
        Arg::new("sigma_deg")
        .short('s')
        .long("sigma")
        .takes_value(true)
        .value_name("sigma in deg")
        .required(true)
        .help("sigma")
    )
    .arg(
        Arg::new("outfile")
        .short('o')
        .long("out")
        .takes_value(true)
        .value_name("output file name")
        .required(true)
        .help("output file name")
    )
    .get_matches();


    let nside=matches.value_of("nside").unwrap().parse::<usize>().unwrap();
    let sigma_deg=matches.value_of("sigma_deg").unwrap().parse::<f64>().unwrap();
    let outfile=matches.value_of("outfile").unwrap();
    let npix=nside2npix(nside);
    let hpmap:Vec<_>=(0..npix)
    .map(|i| if i<npix/2 {(-pix2ang_ring::<f64>(nside, i).pol.to_degrees().powi(2)/(2.0*sigma_deg.powi(2))).exp()}else{0.0}).collect();
    write_map(outfile, &[&hpmap], false, true);
}
