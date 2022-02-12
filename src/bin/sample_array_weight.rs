#![allow(unused_imports)]
extern crate dbf_beam_simulator as lpda;
use clap::{
    App
    , Arg
};

use rand::{
    thread_rng
    , distributions::{
        Uniform
    }
    , Rng
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
    , sym_weight
};

use scorus::{
    mcmc::{
        ensemble_sample::{
            sample
            , init_logprob
            , UpdateFlagSpec
        }
    }
    , linear_space::type_wrapper::LsVec
};

use healpix_fits::{
    read_map
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
    .arg(Arg::new("weights")
        .short('w')
        .long("weights")
        .takes_value(true)
        .value_name("weights")
        .required(true)
        .use_delimiter(true)
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

    //let nants=matches.value_of("nants").unwrap().parse::<usize>().unwrap();
    let init_weights:Vec<_>=matches.values_of("weights").unwrap().map(|x| x.trim().parse::<f64>().unwrap()).collect();
    let nants=init_weights.len()*2+1;
    let spacing=matches.value_of("spacing").unwrap().parse::<f64>().unwrap();
    let target_beam=read_map::<f64>(matches.value_of("target_pattern").unwrap(), &["TEMPERATURE"], 1).pop().unwrap();
    let (averaged_target_beam, _wgt, _theta)=integrate_az(&target_beam);

    let npix=target_beam.len();

    let ant_beam={
        let mut a=read_map::<f64>(matches.value_of("single_antenna_beam").unwrap(), &["TEMPERATURE"], 1);
        a.pop().unwrap()
    };
    assert_eq!(npix, ant_beam.len());
    
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
    let beta=1e6;
    let fobj=|half_w: &LsVec<f64, Vec<f64>>|{
        let weights:Vec<_>=sym_weight(&half_w.0);
        -beta*if weights.iter().any(|&x| x< 0.0 || x>10.0 ){
            f64::INFINITY
        }else{
            let r=beam_opt_func_obj2(&averaged_target_beam, lat, &ant_beam, &ant_x, &ant_y, &ant_z, &weights, &ant_phase, freq);
            //println!("{:?} {:e}", weights, r);
            r
        }
    };

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
    
    let nwalkers=128;

    let mut ensemble:Vec<_>=(0..nwalkers).map(|i|{
        LsVec(if i==0{
            guess.0.clone()
        }else{
            (0..guess.0.len()).map(|_j|{
                rng.sample(Uniform::<f64>::new(0.0, 1.0))
            }).collect()
        })
    }).collect();

    let mut logprob=vec![0.0; nwalkers];
    init_logprob(&fobj, &ensemble, &mut logprob);

    let mut ufs=UpdateFlagSpec::All;
    

    let niter=3000;
    let mut opt_weights=guess.0.clone();
    let mut opt_lp=logprob[0];
    for _i in 0..niter{
        sample(&fobj, &mut ensemble, &mut logprob, &mut rng, 2.0, &mut ufs);
        let optvalue=ensemble.iter().zip(logprob.iter()).reduce(|(e1,lp1),(e2,lp2)|{
            if lp1<lp2{
                (e2,lp2)
            }else{
                (e1,lp1)
            }
        }).unwrap();
        
        if *optvalue.1 > opt_lp{
            opt_lp=*optvalue.1;
            opt_weights=optvalue.0.0.clone();
            println!("{:?} {}", opt_weights, opt_lp/beta);
        }
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
