extern crate dbf_beam_simulator;
use clap::{
    Command, Arg
};

use healpix_fits::{
    read_map
    , write_map
};

use dbf_beam_simulator::{
    utils::integrate_az
    , utils::averaged_beam_to_healpix
};

fn main(){
    let matches=Command::new("mean_by_az")
    .arg(
        Arg::new("input_healpix")
        .short('i')
        .long("input")
        .takes_value(true)
        .value_name("healpix file")
        .required(true)
    )
    .arg(
        Arg::new("output_healpix")
        .short('o')
        .long("output")
        .takes_value(true)
        .value_name("healpix file")
        .required(false)
    )
    .get_matches();

    let hpmap=read_map::<f64>(matches.value_of("input_healpix").unwrap(), &["TEMPERATURE"], 1).pop().unwrap();
    let (mean, wgt, theta)=integrate_az(&hpmap);
    for (&m, (&w, &t)) in mean.iter().zip(wgt.iter().zip(theta.iter())){
        println!("{} {} {}", m, w, t);
    }

    if let Some(fname)=matches.value_of("output_healpix"){
        let hp_data=averaged_beam_to_healpix(&mean);
        write_map(fname, &[&hp_data], false, true);
    }
}
