extern crate dbf_beam_simulator as lpda;
use clap::{
    App
    , Arg
};

use std::{
        io::{
        Write
    }
    , fs::{
        File
    }
};


use serde_yaml::{
    from_reader
    , to_writer
};

use lpda::{
    calc_array_beam
    , calc_phase_from_pointing
    , ArrayCfg
    , AntCfg
};

use healpix_fits::{
    write_map
};

fn main(){
    let matches=App::new("calc_array_beam")
    .arg(Arg::new("array_cfg")
        .short('a')
        .long("array")
        .takes_value(true)
        .value_name("array cfg")
        .required(true)
        .help("array cfg")
    )
    .arg(
        Arg::new("nside")
        .short('s')
        .long("nside")
        .takes_value(true)
        .value_name("nside")
        .required(true)
        .help("nside")
    )
    .arg(
        Arg::new("freq_MHz")
        .short('f')
        .long("freq")
        .takes_value(true)
        .value_name("freq")
        .required(true)
        .help("freq")
    )
    .arg(
        Arg::new("outfile")
        .short('o')
        .long("out")
        .takes_value(true)
        .value_name("outfile")
        .required(true)
        .help("out healpix file")
    )
    .get_matches();


    let mut cfgfile=File::open(matches.value_of("array_cfg").unwrap()).unwrap();
    let nside=matches.value_of("nside").unwrap().parse::<usize>().unwrap();
    let out_file_name=matches.value_of("outfile").unwrap();
    let freq=matches.value_of("freq_MHz").unwrap().parse::<f64>().unwrap()*1e6;
    let array_cfg:ArrayCfg=serde_yaml::from_reader(&mut cfgfile).unwrap();
    let (ant_x, (ant_y, (ant_z, ant_w))):(Vec<f64>, (Vec<f64>, (Vec<f64>, Vec<f64>)))=array_cfg.ants.iter().map(|AntCfg { pos, weight }|{(pos.0,( pos.1,( pos.2, *weight)))}).unzip();
    let ant_phase=calc_phase_from_pointing(&ant_x, &ant_y, &ant_z, 0_f64.to_radians(), 0_f64.to_radians(), 100e6);
    let array_beam=calc_array_beam(nside, &ant_x, &ant_y, &ant_z, &ant_w,&ant_phase, freq, true);
    write_map(out_file_name, &[&array_beam], false, true);
}
