extern crate dbf_beam_simulator as lpda;
use clap::{
    Command
    , Arg
};

use std::{
    fs::{
        File
    }
};


use lpda::{
    calc_phase_from_pointing
    , ArrayCfg
    , AntCfg
    , calc_averaged_array_beam
    , averaged_beam_to_healpix
};


use healpix_fits::{
    write_map
    , read_map
};

fn main(){
    let matches=Command::new("calc_array_beam_total")
    .arg(Arg::new("single_antenna_beam")
        .short('s')
        .long("ant")
        .takes_value(true)
        .value_name("ant healpix")
        .required(true)
        .help("healpix")
    )
    .arg(Arg::new("array_cfg")
        .short('a')
        .long("array")
        .takes_value(true)
        .value_name("array cfg")
        .required(true)
        .help("array cfg")
    )
    .arg(
        Arg::new("freq_MHz")
        .short('f')
        .long("freq")
        .takes_value(true)
        .value_name("freq in MHz")
        .required(true)
        .help("freq in MHz")
    )
    .arg(
        Arg::new("outfile")
        .short('o')
        .long("out")
        .takes_value(true)
        .value_name("outfile")
        .required(false)
        .help("out healpix file")
    )
    .arg(
        Arg::new("lat")
        .short('l')
        .long("lat")
        .takes_value(true)
        .value_name("lat in deg")
        .required(true)
        .help("lat")
    )
    .get_matches();


    let ant_beam={
        let mut a=read_map::<f64>(matches.value_of("single_antenna_beam").unwrap(), &["TEMPERATURE"], 1);
        a.pop().unwrap()
    };
    let mut cfgfile=File::open(matches.value_of("array_cfg").unwrap()).unwrap();
    let lat=matches.value_of("lat").unwrap().parse::<f64>().unwrap();
    
    let freq=matches.value_of("freq_MHz").unwrap().parse::<f64>().unwrap()*1e6;
    let array_cfg:ArrayCfg=serde_yaml::from_reader(&mut cfgfile).unwrap();
    let (ant_x, (ant_y, (ant_z, ant_w))):(Vec<f64>, (Vec<f64>, (Vec<f64>, Vec<f64>)))=array_cfg.ants.iter().map(|AntCfg { pos, weight }|{(pos.0,( pos.1,( pos.2, *weight)))}).unzip();

    let ant_phase=calc_phase_from_pointing(&ant_x, &ant_y, &ant_z, 0_f64.to_radians(), 0_f64.to_radians(), freq);
    
    let (mean, wgt, theta)=calc_averaged_array_beam(lat, &ant_beam, &ant_x, &ant_y, &ant_z, &ant_w, &ant_phase, freq);
    for (&m, (&w, &t)) in mean.iter().zip(wgt.iter().zip(theta.iter())){
        println!("{} {} {}", m, w, t);
    }


    match matches.value_of("outfile"){
        Some(out_file_name)=>{
            let hp_data=averaged_beam_to_healpix(&mean);
            write_map(out_file_name, &[&hp_data], false, true);
        }
        None=>()
    }
}