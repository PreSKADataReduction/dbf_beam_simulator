#![allow(unused_imports)]
extern crate dbf_beam_simulator as lpda;
use clap::{
    Command
    , Arg
};

use rand::{
    thread_rng
};

use std::{
    fs::{
        File
    }
};


use serde_yaml::{
    to_writer
};

use lpda::{
    calc_phase_from_pointing
    , ArrayCfg
    , AntCfg
    , beam_opt_func_obj2
    , integrate_az
    , zenith_ns_sym_array
    , calc_averaged_ant_output
    , calc_averaged_ant_output2
    , calc_array_beam
    , sym_weight
};

use scorus::{
    opt::{
        powell::fmin
        , tolerance::Tolerance
        , pso::{
            Particle
            , ParticleSwarmMaximizer
        }
    }
    , coordinates::{
        rotation3d::RotMatrix
        , Vec3d
    }
    , linear_space::type_wrapper::LsVec
    , healpix::{
        npix2nside
        , rotation::{
            rotate_ring
        }
    }
};

use healpix_fits::{
    read_map
    , write_map
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
    .arg(Arg::new("sky")
        .short('y')
        .long("sky")
        .takes_value(true)
        .value_name("sky healpix")
        .required(true)
        .help("healpix")
    )
    .arg(Arg::new("weights")
        .short('w')
        .long("weights")
        .takes_value(true)
        .value_name("weights")
        .required(true)
        .use_value_delimiter(true)
        .value_delimiter(',')
        .allow_hyphen_values(true)
        .help("number of ants, must be odd")
    )
    .arg(Arg::new("spacing")
        .short('d')
        .long("spacing")
        .takes_value(true)
        .value_name("spacing in meter")
        .required(true)
        .help("antenna spacing")
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
        Arg::new("target_pattern")
        .short('t')
        .long("target")
        .takes_value(true)
        .value_name("healpix")
        .required(true)
        .help("target beam healpix file")
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
    .arg(Arg::new("outfile")
        .short('o')
        .long("out")
        .takes_value(true)
        .value_name("output file name")
        .required(true)
        .help("output file name")
    )
    .get_matches();

    let npart=256;
    //let nants=matches.value_of("nants").unwrap().parse::<usize>().unwrap();
    let init_weights:Vec<_>=matches.values_of("weights").unwrap().map(|x| x.trim().parse::<f64>().unwrap()).collect();
    let nants=init_weights.len()*2+1;
    let spacing=matches.value_of("spacing").unwrap().parse::<f64>().unwrap();
    let target_beam=read_map::<f64>(matches.value_of("target_pattern").unwrap(), &["TEMPERATURE"], 1).pop().unwrap();
    let (averaged_target_beam, _wgt, _theta)=integrate_az(&target_beam);

    let npix=target_beam.len();
    let nside=npix2nside(npix);
    let ant_beam={
        let mut a=read_map::<f64>(matches.value_of("single_antenna_beam").unwrap(), &["TEMPERATURE"], 1);
        a.pop().unwrap()
    };
    assert_eq!(npix, ant_beam.len());
    let sky=read_map::<f64>(matches.value_of("sky").unwrap(), &["TEMPERATURE"], 1).pop().unwrap();
    let target_output=calc_averaged_ant_output(&averaged_target_beam, &sky);
    let lat=matches.value_of("lat").unwrap().parse::<f64>().unwrap();
    
    let freq=matches.value_of("freq_MHz").unwrap().parse::<f64>().unwrap()*1e6;
    
    let y_list:Vec<_>=(1..=(nants-1)/2).map(|i| i as f64*spacing).collect();

    let (ant_x, ant_y, ant_z, ant_w)=zenith_ns_sym_array(&y_list, &init_weights);

    let array_cfg=ArrayCfg{
        ants:ant_x.iter().zip(ant_y.iter().zip(ant_z.iter().zip(ant_w.iter()))).map(|(&x, (&y, (&z, &w)))|{
        AntCfg{
            pos:(x,y,z)
            , weight: w
        }
    }).collect()};

    let mut dump_cfg_file=File::create("array_cfg_dump.yaml").unwrap();
    to_writer(&mut dump_cfg_file, &array_cfg).unwrap();

    let ant_phase=calc_phase_from_pointing(&ant_x, &ant_y, &ant_z, 0_f64.to_radians(), 0_f64.to_radians(), freq);
    
    let diff=beam_opt_func_obj2(&averaged_target_beam, lat, &ant_beam, &ant_x, &ant_y, &ant_z, &ant_w, &ant_phase, freq);
    println!("{}", diff);

    let fobj=|half_w: &LsVec<f64, Vec<f64>>|{
        let weights:Vec<_>=sym_weight(&half_w.0);
        -if weights.iter().any(|&x| x< -10.0 || x>10.0 ){
            f64::INFINITY
        }else{
            let r=beam_opt_func_obj2(&averaged_target_beam, lat, &ant_beam, &ant_x, &ant_y, &ant_z, &weights, &ant_phase, freq);
            //println!("{:?} {:e}", weights, r);
            r
        }
    };

    let rot=RotMatrix::about_axis_by_angle(&Vec3d::new(0.0, 1.0, 0.0), (90.0-lat).to_radians())
    *RotMatrix::about_axis_by_angle(&Vec3d::new(0.0, 0.0, 1.0), -90_f64.to_radians())
    ;

    /*
    let (opt_weights, result)=fmin(&fobj, 
        //&LsVec(vec![1.0; (nants-1)/2])
        &LsVec(init_weights)
        , Tolerance::Abs(1e-10), 1000, Some(&mut |_x, _y|{
        //println!("{:?} {}", x.0, y)
    }));
    */

    let mut rng=thread_rng();
    let guess=LsVec(init_weights.clone());

    println!("init diff={:e}", fobj(&guess));
    
    let mut pso_solver=ParticleSwarmMaximizer::new(&fobj, &LsVec(vec![0.0; (nants-1)/2]),&LsVec(vec![1.0; (nants-1)/2]) , Some(guess), npart, &mut rng);
    let mut opt_weights=init_weights.clone();
    while !pso_solver.converged(0.7, 1e-9, 1e-9) {
        if let Some(ref gbest) = pso_solver.gbest {
            

            opt_weights=gbest.position.0.clone();

            let weights:Vec<_>=sym_weight(&opt_weights);
            let ant_out=calc_averaged_ant_output2(&sky, lat, &ant_beam, &ant_x, &ant_y, &ant_z, &weights, &ant_phase, freq);
            eprintln!("{:?} {:e} {:e}", gbest.position.0, gbest.fitness, (ant_out-target_output)/target_output);
            let array_beam=calc_array_beam(nside, &ant_x, &ant_y, &ant_z, &weights,&ant_phase, freq, true);
            let total_beam:Vec<_>=array_beam.iter().zip(ant_beam.iter()).map(|(&a, &b)|{a*b}).collect();
            let mut rotated_beam=rotate_ring(&total_beam, &rot);
            let max_value=rotated_beam.iter().cloned().reduce(|x,y|{if x>y{x}else{y}}).unwrap();
            rotated_beam.iter_mut().for_each(|x| *x=*x/max_value);

            write_map("total_beam.fits", &[&rotated_beam], false, true);
            
        } else {
            eprint!(".")
        }
        pso_solver.sample(&mut rng, 0.75, 0.5, 1.);
    }


    let ant_w=sym_weight(&opt_weights);
    let array_cfg_opt=ArrayCfg{
        ants:ant_x.iter().zip(ant_y.iter().zip(ant_z.iter().zip(ant_w.iter()))).map(|(&x, (&y, (&z, &w)))|{
        AntCfg{
            pos:(x,y,z)
            , weight: w
        }
    }).collect()};

    let mut dump_cfg_file=File::create(matches.value_of("outfile").unwrap()).unwrap();
    to_writer(&mut dump_cfg_file, &array_cfg_opt).unwrap();
}
