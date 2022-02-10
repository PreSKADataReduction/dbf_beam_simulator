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

use serde_yaml::{
    from_reader
};


use lpda::{
    square_array_beam
    //, calc_phase_from_pointing
    //, zenith_ns_sym_array
    , SquareArrayCfg
};

use scorus::{
    healpix::{
        utils::{
            npix2nside
            , nside2npix
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
    let matches=App::new("calc_square_array_beam")
    .arg(Arg::new("single_antenna_beam")
        .short('s')
        .long("ant")
        .takes_value(true)
        .value_name("ant healpix")
        .required(false)
        .help("healpix")
    )
    .arg(Arg::new("nside")
        .short('n')
        .long("nside")
        .takes_value(true)
        .value_name("nside")
        .required(false)
        .help("nside")
    )
    .arg(Arg::new("array_cfg")
        .short('a')
        .long("array")
        .takes_value(true)
        .value_name("array cfg")
        .required(false)
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
    .arg(
        Arg::new("lon")
        .short('m')
        .long("lon")
        .takes_value(true)
        .value_name("lon in deg")
        .required(false)
        .help("lon")
    )
    .get_matches();

    let (ant_beam, nside)=
    if let Some(fname)=matches.value_of("single_antenna_beam"){
        let mut a=read_map::<f64>(fname, &["TEMPERATURE"], 1);
        let ant_beam=a.pop().unwrap();
        let nside=npix2nside(ant_beam.len());
        (ant_beam, nside)
    }else{
        let nside=matches.value_of("nside").unwrap().parse::<usize>().unwrap();
        let npix=nside2npix(nside);
        let ant_beam=(0..npix).map(|_| 1.0).collect();
        (ant_beam, nside)
    };
    
    let lat=matches.value_of("lat").unwrap().parse::<f64>().unwrap();
    let lon=if let Some(lon)=matches.value_of("lon"){
        lon.parse::<f64>().unwrap()
    }else{
        0.0
    };
    let out_file_name=matches.value_of("outfile").unwrap();
    let freq=matches.value_of("freq_MHz").unwrap().parse::<f64>().unwrap()*1e6;

    let mut cfg_file=File::open(matches.value_of("array_cfg").unwrap()).unwrap();
    let array_cfg:SquareArrayCfg=   from_reader(&mut cfg_file).unwrap();
    

    let array_beam=square_array_beam(nside, &array_cfg, freq, true);
    let total_beam:Vec<_>=array_beam.iter().zip(ant_beam.iter()).map(|(&a, &b)|{a*b}).collect();


    let rot=
    RotMatrix::about_axis_by_angle(&Vec3d::new(0.0, 0.0, 1.0), lon.to_radians())
    *RotMatrix::about_axis_by_angle(&Vec3d::new(0.0, 1.0, 0.0), (90.0-lat).to_radians())
    *RotMatrix::about_axis_by_angle(&Vec3d::new(0.0, 0.0, 1.0), -90_f64.to_radians())
    ;

    let mut rotated_beam=rotate_ring(&total_beam, &rot);

    let max_value=rotated_beam.iter().cloned().reduce(|x,y|{if x>y{x}else{y}}).unwrap();
    rotated_beam.iter_mut().for_each(|x| *x=*x/max_value);

    write_map(out_file_name, &[&rotated_beam], false, true);
}