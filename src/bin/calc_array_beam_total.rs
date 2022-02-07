extern crate dbf_beam_simulator as lpda;
use clap::{
    App
    , Arg
};

use std::{
    fs::{
        File
    }
};



use lpda::{
    calc_array_beam
    , calc_phase_from_pointing
    , ArrayCfg
    , AntCfg
};

use scorus::{
    healpix::{
        utils::{
            npix2nside
        }
        , rotation::rotate_ring
    }
    , coordinates::{
        rotation3d::{
            RotMatrix
        }
        , vec3d::Vec3d
    }
};

use healpix_fits::{
    write_map
    , read_map
};

fn main(){
    let matches=App::new("calc_array_beam_total")
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
        .required(true)
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
    let nside=npix2nside(ant_beam.len());
    let out_file_name=matches.value_of("outfile").unwrap();
    let freq=matches.value_of("freq_MHz").unwrap().parse::<f64>().unwrap()*1e6;
    let array_cfg:ArrayCfg=serde_yaml::from_reader(&mut cfgfile).unwrap();
    let (ant_x, (ant_y, (ant_z, ant_w))):(Vec<f64>, (Vec<f64>, (Vec<f64>, Vec<f64>)))=array_cfg.ants.iter().map(|AntCfg { pos, weight }|{(pos.0,( pos.1,( pos.2, *weight)))}).unzip();

    let ant_phase=calc_phase_from_pointing(&ant_x, &ant_y, &ant_z, 0_f64.to_radians(), 0_f64.to_radians(), freq);
    let array_beam=calc_array_beam(nside, &ant_x, &ant_y, &ant_z, &ant_w,&ant_phase, freq, true);
    let total_beam:Vec<_>=array_beam.iter().zip(ant_beam.iter()).map(|(&a, &b)|{a*b}).collect();


    let rot=RotMatrix::about_axis_by_angle(&Vec3d::new(0.0, 1.0, 0.0), (90.0-lat).to_radians())
    *RotMatrix::about_axis_by_angle(&Vec3d::new(0.0, 0.0, 1.0), -90_f64.to_radians())
    ;

    let rotated_beam=rotate_ring(&total_beam, &rot);

    write_map(out_file_name, &[&rotated_beam], false, true);
}
