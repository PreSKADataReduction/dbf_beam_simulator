extern crate dbf_beam_simulator;

use clap::{
    Command
    ,Arg
    , ArgGroup
};

use rand::{
    thread_rng
};

use ndarray::{
    Ix2
};

use fitsimg::{
    write_img
    , read_img
};


use scorus::{
    healpix::{
        utils::{
            npix2nside
            , nside2npix
        }
    }
    , opt::{
        pso::ParticleSwarmMaximizer
    }
    , linear_space::type_wrapper::LsVec
};

use healpix_fits::{
    read_map
};

use dbf_beam_simulator::{
    quarter_wgt2pattern
    , full2quarter
    , deflattern_quarter_wgt, quarter2full
};



fn main(){
    let matches=Command::new("wgt2beam")
    .arg(
        Arg::new("ant_beam")
        .short('a')
        .long("ant")
        .takes_value(true)
        .value_name("single antenna beam")
        .required(false)
        .help("ant beam")
    )
    .arg(
        Arg::new("target_beam")
        .short('t')
        .long("tb")
        .takes_value(true)
        .value_name("target beam")
        .required(true)
        .help("target_beam")
    )
    .arg(
        Arg::new("wgt0")
        .short('w')
        .long("wgt0")
        .takes_value(true)
        .value_name("wgt file")
        .required(true)
        .help("initial wgt")
    )
    .arg(
        Arg::new("nside")
        .short('n')
        .long("nside")
        .takes_value(true)
        .value_name("nside")
        .required(false)
        .help("nside")
    )
    .arg(
        Arg::new("spacing")
        .short('d')
        .long("spacing")
        .takes_value(true)
        .value_name("spacing in metre")
        .required(true)
        .help("spacing in metre")
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
        Arg::new("sky")
        .short('s')
        .long("sky")
        .takes_value(true)
        .value_name("sky file")
        .required(true)
        .help("sky")
    )
    .arg(
        Arg::new("npart")
        .short('p')
        .long("npart")
        .takes_value(true)
        .value_name("num of particles")
        .required(false)
        .help("num of particles")
    )
    .arg(
        Arg::new("out_wgt")
        .short('o')
        .long("out")
        .takes_value(true)
        .value_name("out_wgt")
        .required(true)
        .help("out_wgt")
    )
    .group(
        ArgGroup::new("ant_beam_or_nside")
        .args(&["ant_beam", "nside"])
        .required(true)
    )
    .get_matches();

    let (ant_beam, nside)=
    if let Some(fname)=matches.value_of("ant_beam"){
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

    let npix=nside2npix(nside);
    let sky=read_map::<f64>(matches.value_of("sky").unwrap(), &["TEMPERATURE"], 1).pop().unwrap();
    assert_eq!(npix, sky.len());
    
    let target_beam=read_map::<f64>(
        matches.value_of("target_beam").unwrap(), &["TEMPERATURE"], 1).pop().unwrap();
    let norm=target_beam.iter().cloned().sum::<f64>();
    let target_ant_out=target_beam.iter().zip(sky.iter()).map(|(&b,&s)| b*s).sum::<f64>()/norm;

    let target_beam:Vec<_>=target_beam.iter().map(|&x| x/norm).collect();
    assert_eq!(npix, target_beam.len());

    
    
    
    let wgt=read_img::<f64>(matches.value_of("wgt0").unwrap().to_string(), 0).unwrap().into_dimensionality::<Ix2>().unwrap();
    let h=wgt.shape()[0];
    let w=wgt.shape()[1];

    let d=matches.value_of("spacing").unwrap().parse::<f64>().unwrap();

    let freq_mhz=matches.value_of("freq_MHz").unwrap().parse::<f64>().unwrap();

    let wgt_eff=full2quarter(wgt.view());
    

    
    let fobj=|x: &LsVec<f64, Vec<f64>>|{
        let wgt=deflattern_quarter_wgt(&x.0, h, w);
        let array_beam=quarter_wgt2pattern(
            wgt.view()
            , d, freq_mhz, nside);
        let total_beam:Vec<_>=array_beam.iter().zip(ant_beam.iter()).map(|(&a,&b)| a*b).collect();
        let norm=total_beam.iter().sum::<f64>();
        -(total_beam.iter().zip(target_beam.iter()).map(|(&x,&y)| (x/norm-y).powi(2)).sum::<f64>()*npix as f64).log10()
    };

    let guess=LsVec::<f64, Vec<_>>(wgt_eff.iter().skip(1).cloned().collect());
    let npart=if let Some(np)=matches.value_of("npart"){
        np.parse().unwrap()
    }else{
        64
    };
    let mut rng=thread_rng();

    eprintln!("init diff::{}", fobj(&guess));
    let mut opt_weights=guess.0.clone();
    
    let mut pso_solver=
    ParticleSwarmMaximizer::new(
        &fobj, &LsVec(vec![0.0; h*w/4-1]),&LsVec(vec![1.0; h*w/4-1]) ,
        Some(guess), npart, &mut rng);
    
    while !pso_solver.converged(0.7, 1e-9, 1e-9) {
        if let Some(ref gbest) = pso_solver.gbest {
            opt_weights=gbest.position.0.clone();

            let wgt=deflattern_quarter_wgt(&opt_weights, h, w);
            let array_beam=quarter_wgt2pattern(
                wgt.view()
                , d, freq_mhz, nside);
            let total_beam:Vec<_>=array_beam.iter().zip(ant_beam.iter()).map(|(&a,&b)| a*b).collect();
            let norm=total_beam.iter().sum::<f64>();
            let ant_out=total_beam.iter().zip(sky.iter()).map(|(&b,&s)| b*s).sum::<f64>()/norm;
            eprintln!("{} {} {}", gbest.fitness,target_ant_out, (ant_out-target_ant_out).abs()/target_ant_out);
        } else {
            eprint!(".")
        }
        pso_solver.sample(&mut rng, 0.75, 0.5, 1.);
    }

    let wgt=quarter2full(deflattern_quarter_wgt(&opt_weights, h, w).view());
    write_img(matches.value_of("out_wgt").unwrap().to_string(), &wgt.into_dyn()).unwrap();
}
