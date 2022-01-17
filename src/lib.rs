use std::{
    f64::{
        consts::{
            PI
        }
    }
};

use serde::{Serialize, Deserialize};


use num::{
    complex::{
        Complex
    }
};

const LIGHT_SPEED:f64=2.99792458E8;

use scorus::{
    healpix::{
        utils::{
            nside2npix
            , npix2nside
            , nside2nring
        }
        , pix::{
            pix2vec_ring, pix2ring_ring
            , ring2z_ring
        }
    }
    , coordinates::{
        SphCoord
        ,Vec3d
    }
};

#[derive(Clone, Serialize, Deserialize)]
pub struct ArrayCfg{
    pub ants: Vec<AntCfg>
}

#[derive(Clone, Serialize, Deserialize)]
pub struct AntCfg{
    pub pos: (f64,f64,f64)
    , pub weight: f64
}

pub fn calc_array_beam(nside: usize, x_list: &[f64], y_list: &[f64], z_list: &[f64], w_list: &[f64], phi_list:&[f64], freq_Hz: f64, ground_cut: bool)->Vec<f64>{
    let npix=nside2npix(nside);
    let lambda=LIGHT_SPEED/(freq_Hz);
    (0..npix).map(|i|{
        if i<npix/2 || !ground_cut{
            let pointing=pix2vec_ring::<f64>(nside, i);
            x_list.iter().zip(y_list.iter().zip(z_list.iter().zip(w_list.iter().zip(phi_list.iter())))).map(|(&x, (&y, (&z, (&w, &phi)))) |{
                let dl=pointing[0]*x+pointing[1]*y+pointing[2]*z;
                let phase=dl/lambda*2.0*PI;

                Complex::from_polar(w, phase-phi)
            }).sum::<Complex<f64>>().norm_sqr()
        }else{
            0.0
        }
    }).collect()
}

pub fn calc_phase_from_pointing(x_list: &[f64], y_list: &[f64], z_list: &[f64], az_from_north: f64, zenith: f64, freq_Hz: f64)->Vec<f64>{
    let az_from_x=PI/2.0-az_from_north;
    let lambda=LIGHT_SPEED/freq_Hz;
    let dir=Vec3d::from_sph_coord(SphCoord::new(zenith, az_from_x));
    x_list.iter().zip(y_list.iter().zip(z_list.iter())).map(|(&x, (&y, &z))|{
        let dl=dir[0]*x+dir[1]*y+dir[2]*z;
        dl/lambda*2.0*PI
    }).collect()
}

pub fn integrate_az(hmap: &[f64])->(Vec<f64>, Vec<usize>, Vec<f64>){
    let npix=hmap.len();
    let nside=npix2nside(npix);
    let nring=nside2nring(nside);
    let mut wgt=vec![0_usize; nring];
    let mut mean_values=vec![0.0; nring];
    let ring_idx:Vec<_>=(0..npix).map(|ipix| pix2ring_ring(nside, ipix)).collect();
    ring_idx.iter().for_each(|&i| wgt[i-1]+=1);
    hmap.iter().zip(ring_idx.iter()).for_each(|(&x, &iring)|{
        mean_values[iring-1]+=x;
    });
    mean_values.iter_mut().zip(wgt.iter()).for_each(|(x, &w)|{
        *x/=w as f64
    });

    let theta=(1..=nring).map(|iring| ring2z_ring::<f64>(nside, iring).acos()).collect();
    (mean_values, wgt, theta)
}
