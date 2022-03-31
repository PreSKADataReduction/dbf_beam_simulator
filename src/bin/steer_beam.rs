extern crate dbf_beam_simulator as lpda;
use clap::{
    Command
    , Arg
    , ArgGroup
};


use scorus::{
    healpix::{
        rotation::rotate_ring
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
    let matches=Command::new("steer_beam")
    .arg(Arg::new("beam")
        .short('b')
        .long("beam")
        .takes_value(true)
        .value_name("beam healpix")
        .required(false)
        .help("healpix")
    )
    .arg(Arg::new("sky")
        .short('s')
        .long("sky")
        .takes_value(true)
        .value_name("sky")
        .required(false)
        .help("sky healpix")
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
        .allow_hyphen_values(true)
        .help("lat")
    )
    .arg(
        Arg::new("lon")
        .short('m')
        .long("lon")
        .takes_value(true)
        .value_name("lon in deg")
        .required(false)
        .allow_hyphen_values(true)
        .help("lon")
    )
    .group(
        ArgGroup::new("inputs")
        .args(&["beam", "sky"])
        .required(true)
    )
    .get_matches();

    let lat=matches.value_of("lat").unwrap().parse::<f64>().unwrap();
    let lon=if let Some(lon)=matches.value_of("lon"){
        lon.parse::<f64>().unwrap()
    }else{
        0.0
    };

    

    let (hp_data, rot)=
    if let Some(fname)=matches.value_of("beam"){
        let data=read_map::<f64>(fname, &["TEMPERATURE"], 1).pop().unwrap();
        let rot=
        RotMatrix::about_axis_by_angle(&Vec3d::new(0.0, 0.0, 1.0), (90.0+lon).to_radians())
        *RotMatrix::about_axis_by_angle(&Vec3d::new(1.0, 0.0, 0.0), (90.0-lat).to_radians());
        (data, rot)
    }else if let Some(fname)=matches.value_of("sky"){
        let data=read_map::<f64>(fname, &["TEMPERATURE"], 1).pop().unwrap();
        let rot=
        RotMatrix::about_axis_by_angle(&Vec3d::new(1.0, 0.0, 0.0), -(90.0-lat).to_radians())
        *RotMatrix::about_axis_by_angle(&Vec3d::new(0.0, 0.0, 1.0), -(90.0+lon).to_radians());
        (data,rot)
    }else{
        panic!()
    };
    
    println!("{:?}", rot);
    let out_file_name=matches.value_of("outfile").unwrap();
    
    let rotated=rotate_ring(&hp_data, &rot);

    write_map(out_file_name, &[&rotated], false, true);
}
