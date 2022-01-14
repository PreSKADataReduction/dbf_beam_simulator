use std::fs::{
    read_to_string
};

use pest::Parser;

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
    let mut context=parse_nec_file(NecParser::parse(Rule::NecFile, &read_to_string("./21cmax3.nec").unwrap()).unwrap().next().unwrap());

    context.nec_fr_card(0, 1, 120.0, 0.0);
    


    let nside=256;
    let npix=nside2npix(nside);
    let angular_resolution=(4.0*f64::PI()/npix as f64).sqrt().to_degrees();
    println!("{}", angular_resolution);
    let nphi=(360.0/angular_resolution*4.0) as usize+1;
    let ntheta=(nphi-1)/2+1;
    let dphi=360.0/(nphi-1) as f64;
    let dtheta=180.0/(ntheta-1) as f64;
    
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

    write_map("a.fits", &[&data], false, true);
}
