use std::fs::{
    read_to_string
};

use pest::Parser;

use clap::{
    Command
    , Arg
};

use num::{
    traits::{
        FloatConst
    }
};

use healpix_fits::write_map;

use scorus::{
    healpix::{
        interp::{
            get_interpol_ring
        }
        , utils::{
            nside2npix
        }
    }
    , coordinates::{
        SphCoord
    }
};

use necrs::nec_parser::{
    parse_nec_file
    , NecParser
    , Rule
};

pub fn main(){
    let matches=Command::new("calc_ant_beam")
    .arg(Arg::new("nec")
        .short('n')
        .long("nec")
        .takes_value(true)
        .value_name("nec file")
        .required(true)
        .help("nec file")
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

    let nec_file_name=matches.value_of("nec").unwrap();
    let freq=matches.value_of("freq_MHz").unwrap().parse::<f64>().unwrap();
    let nside=matches.value_of("nside").unwrap().parse::<usize>().unwrap();
    let out_file_name=matches.value_of("outfile").unwrap();

    let mut context=parse_nec_file(NecParser::parse(Rule::NecFile, &read_to_string(nec_file_name).unwrap()).unwrap().next().unwrap());

    context.nec_fr_card(0, 1, freq, 0.0);
    

    let npix=nside2npix(nside);
    let angular_resolution=(4.0*f64::PI()/npix as f64).sqrt().to_degrees();
    println!("{}", angular_resolution);
    
    
    
    //context.nec_rp_card(0, ntheta as i32, nphi as i32, 1, 0, 0, 0, 0.0, 0.0, dtheta, dphi, 0.0, 0.0);
    let (thetas, phis)=context.rp_from_npix(npix*4, 0, 1, 0, 0, 0, 0.0, 0.0);
    let mut data=vec![0.0; npix];
    let mut wgt=vec![0.0; npix];

    for (i, &theta) in thetas.iter().enumerate(){
        if theta > 90.0{
            continue
        }
        for (j, &phi) in phis.iter().enumerate(){
            let g=(context.nec_gain(0, i as i32, j as i32)/10.0).exp();
            let dir=SphCoord::new(theta.to_radians(), phi.to_radians());
            let (pix, w)=get_interpol_ring(nside, dir);
            for (&p,&w) in pix.iter().zip(w.iter()){
                wgt[p]+=w;
                data[p]+=w*g;
            }
        }
    }    for (d,&w) in data.iter_mut().zip(wgt.iter()){
        if w>0.0{
            *d/=w;
        }
    }

    write_map(out_file_name, &[&data], false, true);
}
