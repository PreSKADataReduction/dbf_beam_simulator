use scorus::{
    healpix::{
        pix::{
            pix2ring_ring
            , ring2z_ring
        }
        , utils::{
            nside2npix
        }
    }
};

use healpix_fits::{
    write_map
};

extern crate lpda_beam_healpix;

use lpda_beam_healpix::{
    integrate_az
};


fn main(){
    let nside=16;
    let npix=nside2npix(nside);
    let map:Vec<_>=(0..npix).map(|ipix| 1.0 ).collect();

    let (m, w, theta)=integrate_az(&map);

    for (&m1, (&w1, &theta1)) in m.iter().zip(w.iter().zip(theta.iter())){
        println!("{} {} {}", m1, w1, theta1);
    }
}