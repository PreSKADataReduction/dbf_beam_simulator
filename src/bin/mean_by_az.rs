extern crate dbf_beam_simulator;
use clap::{
    App, Arg
};

use healpix_fits::{
    read_map
};

use dbf_beam_simulator::{
    integrate_az
};

fn main(){
    let matches=App::new("mean_by_az")
    .arg(
        Arg::new("input_healpix")
        .short('i')
        .long("input")
        .takes_value(true)
        .value_name("healpix file")
        .required(true)
    ).get_matches();

    let hpmap=read_map::<f64>(matches.value_of("input_healpix").unwrap(), &["TEMPERATURE"], 1).pop().unwrap();
    let (mean, wgt, theta)=integrate_az(&hpmap);
    for (&m, (&w, &t)) in mean.iter().zip(wgt.iter().zip(theta.iter())){
        println!("{} {} {}", m, w, t);
    }
}